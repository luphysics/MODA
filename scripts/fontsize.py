"""
Python script which modifies the font sizes of all MODA's plots directly in the code.

This is the origin of these lines of code in MODA source files: 
    globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.

To modify the font sizes in MODA to a new value X, open a terinal in the MODA folder and run:
    python scripts/fontsize.py X

For example, for a font size of 12:
    python scripts/fontsize.py 12

Note: On macOS/Linux you may need to use `python3` instead of `python`.
"""

import os
from os import path
from pathlib import Path
import sys
from typing import List

assert sys.version_info >= (3, 6), "Please run with Python 3.6 or above."

args = sys.argv[1:]
assert (
    args
), "Font size must be supplied as an argument. For example, 'python fontsize.py 1' for a font size of 12."

# Set working directory to MODA folder.
here = path.abspath(path.dirname(__file__))
os.chdir(Path(here).parent)

gfs = "globalfontsize"
comment = "Do not edit this line manually. See scripts/fontsize.py."


def readlines(filename: str) -> List[str]:
    with open(filename, "r") as f:
        lines = f.readlines()

    return lines


def mutate_lines(filename: str, lines: List[str]) -> List[str]:
    lines = lines.copy()
    for index, item in enumerate(lines.copy()):
        if item.strip().startswith(gfs):
            indent = get_indentation(lines[index + 1 :])
            new = f"{indent}{gfs} = {fontsize}; % {comment}\n"

            lines[index] = new
            print(f"{filename}: editing line {index+1}.")

    return lines


def get_indentation(lines: List[str]) -> str:
    for l in lines:
        if l.strip():
            return " " * (len(l) - len(l.lstrip()))
    return ""


def write_lines(filename: str, lines: List[str]) -> None:
    with open(filename, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    fontsize = int(args[0])

    for filepath in Path(".").rglob("*.m"):
        filename = filepath.absolute()

        lines = readlines(filename)
        mutated = mutate_lines(path.split(filename)[1], lines)

        if lines != mutated:
            write_lines(filename, mutated)

    print(f"Fonts updated to size {fontsize}.")
