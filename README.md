# computational-group-theory
Project in Computational Group Theory for COMPSCI 380 at The University of Auckland.

## Supervisors
[Cris Calude](https://www.cs.auckland.ac.nz/~cristian/) and [Eamonn O'Brien](https://www.math.auckland.ac.nz/~obrien/)

## Project Goals and Outline
The rough project outline is (currently) as follows:

  1. Complete a significant number of exercises of the book Presentations of Groups by D. L. Johnson.
  COMPLETED. The majority of the exercises from the first four chapters have been completed. I have read further than this of course. In the future I may come back and complete more of the exercises.

  2. Write an implementation of an algorithm of the computation of the Smith Normal Form of an matrix with integer entries. Such a procedure was discussed in Maths 720, as it can be used to compute the factors described in the fundamental theorem of finitely generated abelian groups. COMPLETED. See `magma` folder and report.

  3. Read relevant literature on implementation-specific considerations around the algorithm implemented in (2), for example, how to deal with integer overflow. Write a report on these considerations and improve the implementation using the knowledge gained from this. Run the algorithm on non-trivial test cases. COMPLETED. See report.

This project has been submitted for assessment as at 30th October 2017.

## Usage
Clone this repository. To read the solutions to the exercises, compile `exercises/johnson.tex` using your preferred LaTeX distribution. To read the report, compile `report/report.tex`. All of the magma code is included in the folder `magma`.

