"""Verified, resumable-by-archive retrieval and strict pairing of public images.

Only verified complete files are reused. Raw ZIPs are not extracted into arbitrary
paths; selected image/mask members are read directly into memory by the worker.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import ssl
import time
from urllib.request import Request, urlopen
from urllib.error import URLError
import zipfile
IMAGE_DIRS = {'images': 'train', 'val_images': 'val', 'test_images': 'test'}
MASK_DIRS = {'masks': 'train', 'val_masks': 'val', 'test_masks': 'test'}
EXTENSIONS = {'.tif', '.tiff', '.png', '.bmp', '.jpg', '.jpeg'}

def hash_file(path: Path, algorithm='sha256') -> str:
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1 << 20), b''):
            h.update(data)
    return h.hexdigest()

def validate_member(name: str) -> None:
    p = PurePosixPath(name)
    if '\\' in name or p.is_absolute() or '..' in p.parts or any((':' in s for s in p.parts)) or ('\x00' in name):
        raise ValueError('Unsafe archive member: ' + repr(name))

def verify_archive(path: Path, md5: str) -> dict:
    if hash_file(path, 'md5') != md5:
        raise ValueError('Archive MD5 differs from the pinned publisher checksum: ' + path.name)
    with zipfile.ZipFile(path) as z:
        total = 0
        seen = set()
        for info in z.infolist():
            validate_member(info.filename)
            if info.filename in seen:
                raise ValueError('Duplicate ZIP member')
            seen.add(info.filename)
            total += info.file_size
            if info.external_attr >> 16 & 61440 == 40960:
                raise ValueError('Symlink in data archive')
            if info.file_size > 500000000:
                raise ValueError('Unexpectedly large individual archive member')
        if total > 2000000000:
            raise ValueError('Data archive expands beyond the declared safety limit')
        bad = z.testzip()
        if bad:
            raise ValueError('ZIP CRC failed: ' + bad)
    return {'file': path.name, 'md5': md5, 'sha256': hash_file(path), 'bytes': path.stat().st_size, 'archive_members': len(seen), 'expanded_bytes': total}

def retrieve(collection: dict, record_id: str, cache: Path, downloads: Path | None=None, timeout=90) -> tuple[Path, dict]:
    cache.mkdir(parents=True, exist_ok=True)
    destination = cache / collection['file']
    md5 = collection['md5']
    errors = []
    if destination.exists():
        try:
            receipt = verify_archive(destination, md5)
            receipt.update(acquisition='verified_cache')
            return (destination, receipt)
        except Exception as e:
            errors.append('cached file not reused: ' + str(e))
            destination.rename(destination.with_name(destination.name + '.invalid_' + str(time.time_ns())))
    if downloads and downloads.exists():
        pattern = re.compile('^' + re.escape(Path(collection['file']).stem) + '(?:\\s*\\(\\d+\\))?\\.zip$', re.I)
        for candidate in sorted((p for p in downloads.iterdir() if p.is_file() and pattern.match(p.name)), key=lambda x: x.stat().st_mtime, reverse=True):
            if hash_file(candidate, 'md5') == md5:
                shutil.copy2(candidate, destination)
                receipt = verify_archive(destination, md5)
                receipt.update(acquisition='verified_downloads_copy', source_filename=candidate.name)
                return (destination, receipt)
    urls = [f"https://zenodo.org/records/{record_id}/files/{collection['file']}?download=1", f"https://zenodo.org/api/records/{record_id}/files/{collection['file']}/content"]
    temp = destination.with_suffix('.partial')
    for url in urls:
        for attempt in range(2):
            try:
                print('DOWNLOAD: ' + collection['file'] + '; attempt ' + str(attempt + 1), flush=True)
                request = Request(url, headers={'User-Agent': 'RAPD-source-transfer/1.0 (scientific-data-retrieval)', 'Accept': 'application/octet-stream'})
                with urlopen(request, timeout=timeout) as response, temp.open('wb') as stream:
                    received = 0
                    last = 0
                    start = time.monotonic()
                    content_type = response.headers.get('Content-Type', '')
                    if 'text/html' in content_type:
                        raise ValueError('Server returned HTML instead of a ZIP')
                    while True:
                        block = response.read(1 << 20)
                        if not block:
                            break
                        stream.write(block)
                        received += len(block)
                        if received > 500000000:
                            raise ValueError('Unexpected download size')
                        if time.monotonic() - start > 1800:
                            raise TimeoutError('Download exceeded 30 minutes')
                        if received - last >= 10000000:
                            print(f'  downloaded {received / 1000000.0:.1f} MB', flush=True)
                            last = received
                receipt = verify_archive(temp, md5)
                temp.replace(destination)
                receipt.update(file=destination.name, acquisition='https_download', url=url, prior_attempt_errors=errors)
                return (destination, receipt)
            except Exception as e:
                errors.append(type(e).__name__ + ': ' + str(e))
                print('DOWNLOAD_RETRY: ' + errors[-1], flush=True)
                if temp.exists():
                    temp.unlink()
                if attempt == 0:
                    time.sleep(2)
    raise RuntimeError('Cannot retrieve verified ' + collection['file'] + '. Details: ' + ' | '.join(errors))

def pair_members(path: Path) -> tuple[list[dict], list[dict]]:
    images = {}
    masks = {}
    notes = []
    with zipfile.ZipFile(path) as z:
        for info in z.infolist():
            name = info.filename
            validate_member(name)
            if info.is_dir():
                continue
            p = PurePosixPath(name)
            if '__MACOSX' in p.parts or any((s.startswith('.') for s in p.parts)):
                continue
            if p.suffix.lower() not in EXTENSIONS:
                continue
            hits = [(i, part.lower()) for i, part in enumerate(p.parts[:-1]) if part.lower() in IMAGE_DIRS | MASK_DIRS]
            if len(hits) != 1:
                notes.append({'member': name, 'status': 'UNRECOGNIZED_IMAGE_OR_MASK_DIRECTORY'})
                continue
            i, folder = hits[0]
            table = images if folder in IMAGE_DIRS else masks
            partition = (IMAGE_DIRS | MASK_DIRS)[folder]
            parent = '/'.join(p.parts[:i])
            rel = PurePosixPath(*p.parts[i + 1:]).with_suffix('').as_posix()
            key = (parent.casefold(), partition, rel.casefold())
            if key in table:
                raise ValueError('Ambiguous duplicate image/mask pairing: ' + name)
            table[key] = name
        pairs = []
        for key in sorted(images):
            if key not in masks:
                notes.append({'member': images[key], 'status': 'IMAGE_WITHOUT_MASK'})
                continue
            image_name = images[key]
            mask_name = masks[key]
            pairs.append({'image_member': image_name, 'mask_member': mask_name, 'source_partition': key[1], 'id': hashlib.sha256((path.stem + '\n' + image_name).encode()).hexdigest()[:16], 'original_name': PurePosixPath(image_name).name})
        for key in sorted(set(masks) - set(images)):
            notes.append({'member': masks[key], 'status': 'MASK_WITHOUT_IMAGE'})
    return (pairs, notes)

def try_metadata(record_id: str, output: Path) -> dict:
    """Metadata failure is recorded, not confused with a failed verified data download."""
    url = f'https://zenodo.org/api/records/{record_id}'
    try:
        with urlopen(Request(url, headers={'User-Agent': 'RAPD-source-transfer/1.0'}), timeout=20) as r:
            data = r.read(5000001)
        if len(data) > 5000000:
            raise ValueError('Unexpected metadata size')
        obj = json.loads(data)
        if str(obj.get('id')) != str(record_id):
            raise ValueError('Record identifier mismatch')
        output.write_bytes(data)
        return {'status': 'RETRIEVED', 'url': url, 'license_metadata': obj.get('metadata', {}).get('license'), 'sha256': hash_file(output)}
    except Exception as e:
        return {'status': 'METADATA_NOT_RETRIEVED', 'url': url, 'error': str(e), 'archive_checksums_remain_pinned': True}
