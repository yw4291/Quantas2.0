import sys
import json

def parse_config(config_file):
    groups = {}
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            group, sample = line.split()
            if group not in groups:
                groups[group] = []
            groups[group].append(sample)

    with open("samples.json", "w") as out:
        json.dump(groups, out, indent=4)

if __name__ == "__main__":
    parse_config(sys.argv[1])
