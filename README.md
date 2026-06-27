# My C Protein Interpreter
A simple (unsafe) program that reads a .faa file and creates a .csv with each protein name, and various specified properties


**You likely shouldn't use this; it was made over the span of a couple of days, and is neither safe nor effective.**

## Example usage 
> program.exe -f [input .faa file] -o [output .csv file] -mw (if you want the molecular weight) -im (if you want a guess if the protein is integral) -nc [pH] (if you want the net charge at a given pH)

 * the input file must be have this layout

> > [Identifier][SPACE].....
> [AA Sequence]

## How to build

> make -f Makefile

* be sure to specify the compiler version  and the compile location you want in the makefile

