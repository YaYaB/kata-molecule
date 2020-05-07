import json
import time
from collections import defaultdict

# Some parameters
beg_separators = ["(", "[", "{"]
end_separators = {"(": ")", "[": "]", "{": "}"}


def update(prev, prev_number, tmp_dict, output):
    if prev == "":
        return
    nb = 1
    if prev_number:
        nb = int(prev_number) 
    # If it was not a submolecule
    if not tmp_dict:
        output[prev] += nb
    else: # If it was a submolecule
        for elem, val in tmp_dict.items():
            output[elem] += val * nb


def parse_molecule(molecule, ending=None):
    output = defaultdict(lambda: 0)
    prev, prev_number, tmp_dict = "", "", None
    i = 0
    len_molecule = len(molecule)
    while i < len_molecule:
        c = molecule[i]
        #print(c)
        # Check if it is an uppercase
        if c.isupper():
            update(prev, prev_number, tmp_dict, output)
            prev, prev_number, tmp_dict = c, "", None
        # Check if it s a lowercase
        elif c.islower():
            prev += c
        # Check if it is a digit
        elif c.isdigit():
            prev_number += c
        # Check if it is a beggining separator
        elif c in beg_separators:
            update(prev, prev_number, tmp_dict, output)
            tmp_dict, offset = parse_molecule(molecule[i+1:], end_separators[c])
            if offset == -1:
                raise ValueError('Formula is not correct due to absence of {}'.format(end_separators[c]))
            prev, prev_number = c, ""
            i += offset 
            continue
        
        # Check if it is the ending separator
        elif c == ending:
            update(prev, prev_number, tmp_dict, output)
            return output, i + 2
        # Otherise through an error
        else:
            raise ValueError('Formula is not correct due to {} is not valid'.format(c))
        i += 1
        
    # Final character
    update(prev, prev_number, tmp_dict, output)

    return output, -1


if __name__ == "__main__":

    examples = ['H20', 'FC(Br)(Cl)F', 'CH4', 'N20']

    for molecule in examples:
        print(molecule, dict(parse_molecule(molecule)[0]))
