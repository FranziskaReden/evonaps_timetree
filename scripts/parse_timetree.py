import argparse
from utils import read_timetree

def write_out_tax_names(names, file="data/TimeTree5_leaves.txt"): 

    with open(file, "w") as w: 
        for name in names: 
            w.write(name+"\n")

def main():

    parser = argparse.ArgumentParser(description="Parse a TimeTree file and extract leave information.")
    parser.add_argument("--tree", type=str, help="Path to the tree file to parse.")
    args = parser.parse_args()

    # Read the tree file and extract names
    tree, leaf_names = read_timetree(args.tree)

    # Write out name of the leaves in a txt file
    write_out_tax_names([clade.names for clade in leaf_names], file="data/TimeTree5_leaves.txt")

    return 0

if __name__ == "__main__":
    main()
