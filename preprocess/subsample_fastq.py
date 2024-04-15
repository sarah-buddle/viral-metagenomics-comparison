import argparse

def subsample_fastq(input_fastq, output_fastq, target_bases):
    total_bases = 0
    line_number = 0

    with open(input_fastq, "r") as infile, open(output_fastq, "w") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break

            sequence = infile.readline()
            if not sequence:
                break

            quality_header = infile.readline()
            if not quality_header:
                break

            quality = infile.readline()
            if not quality:
                break

            line_number += 4

            total_bases += len(sequence.strip())
            if total_bases >= target_bases:
                break

            outfile.write(header)
            outfile.write(sequence)
            outfile.write(quality_header)
            outfile.write(quality)

    print(f"Total bases written: {total_bases}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subsample a FASTQ file to a target base count.")
    parser.add_argument("input_fastq", help="Input FASTQ file")
    parser.add_argument("output_fastq", help="Output FASTQ file")
    parser.add_argument("target_bases", type=int, help="Target base count")

    args = parser.parse_args()
    subsample_fastq(args.input_fastq, args.output_fastq, args.target_bases)
