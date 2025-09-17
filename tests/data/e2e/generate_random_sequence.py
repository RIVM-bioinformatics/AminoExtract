"""Quick script to generate a random DNA sequence of a given length."""

import secrets


def generate_random_sequence(length: int, no_stop_codons=False) -> str:
    """Generate a random DNA sequence of a given length."""

    seq = "".join(secrets.choice("ATGC") for _ in range(length))

    if not no_stop_codons:
        return seq

    stop_codons = ["TAA", "TAG", "TGA"]

    def _detector(seq):
        for stop_codon in stop_codons:
            if stop_codon in seq:
                return True
        return False

    def _replacer(seq):
        for stop_codon in stop_codons:
            seq = seq.replace(stop_codon, "".join(secrets.choice("ATGC") for _ in range(3)))
        return seq

    while _detector(seq):
        seq = _replacer(seq)

    return seq


if __name__ == "__main__":
    print(generate_random_sequence(10_000, True))
