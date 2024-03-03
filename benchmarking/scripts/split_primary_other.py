import sys
import os

def main(input_fasta, outdir):
    os.makedirs(outdir, exist_ok=True)
    bname = os.path.splitext(os.path.basename(input_fasta))[0]
    outprimary = os.path.join(outdir, bname + "_primary.fa")
    outother = os.path.join(outdir, bname + "_other.fa")
    with open(input_fasta, "r") as inf, open(outprimary, "w") as outp, open(outother, "w") as outo:
        primary = True
        for i, line in enumerate(inf):
            if i != 0 and line.startswith(">"):
                primary = False
            if primary:
                outp.write(line)
            else:
                outo.write(line)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 2:
        raise ValueError ("USAGE: split_primary_other.py fasta path/to/outdir")
    if os.path.exists(args[1]):
        print("warning: output directory already exists")
    main(args[0], args[1])
