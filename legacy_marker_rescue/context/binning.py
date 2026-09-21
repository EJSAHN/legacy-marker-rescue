"""The same size-clustering rule used in the completed validation run."""

def cluster_sizes(sizes, tolerance):
    """Archived greedy bin rule based on unique integer sizes (not counts)."""
    bins = []
    total = 0
    for size in sorted(set((int(s) for s in sizes))):
        center = total / len(bins[-1]) if bins else 0
        if bins and abs(size - center) / max(center, 1) <= tolerance:
            bins[-1].append(size)
            total += size
        else:
            bins.append([size])
            total = size
    mapping = {}
    centers = []
    for i, members in enumerate(bins):
        centers.append(float(round(sum(members) / len(members))))
        for s in members:
            mapping[s] = i
    return (mapping, centers, bins)
