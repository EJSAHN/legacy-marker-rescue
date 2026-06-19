from legacy_marker_rescue.sequence import reverse_complement, gc_fraction


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGC") == "GCTT"


def test_gc_fraction():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("AATT") == 0.0
