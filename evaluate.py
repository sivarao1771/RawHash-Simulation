import sys

def parse_true(read_name):
    parts = read_name.split('!')
    if len(parts) < 5:
        return None
    ref = parts[1]
    start = int(parts[2])
    end = int(parts[3])
    strand = parts[4]
    return ref, start, end, strand

def overlap(s1, e1, s2, e2):
    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    if overlap_end <= overlap_start:
        return False
    overlap_len = overlap_end - overlap_start
    # Use mapped region length as denominator, not true region
    mapped_len = e2 - s2
    return (overlap_len / mapped_len) >= 0.5

tp, fp, fn = 0, 0, 0

with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        read_name = cols[0]
        mapped = cols[4] != '*'

        true = parse_true(read_name)
        if true is None:
            continue
        true_ref, true_start, true_end, true_strand = true

        if not mapped:
            fn += 1
            continue

        mapped_ref = cols[5]
        mapped_start = int(cols[7])
        mapped_end = int(cols[8])

        correct_ref = (mapped_ref == true_ref)
        correct_overlap = overlap(true_start, true_end, mapped_start, mapped_end)

        if correct_ref and correct_overlap:
            tp += 1
        else:
            fp += 1
            fn += 1

precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall    = tp / (tp + fn) if (tp + fn) > 0 else 0
f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

print(f"TP: {tp}  FP: {fp}  FN: {fn}")
print(f"Precision : {precision:.4f}")
print(f"Recall    : {recall:.4f}")
print(f"F1 Score  : {f1:.4f}")
