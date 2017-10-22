procedure Diagonalise(~A)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    if not IsZero(A) and m ne 0 and n ne 0 then
	// find a nonzero entry
	i := 1; j := 1;
	while A[i,j] eq 0 do
	    if j lt n then
		j +:= 1;
	    else
		i +:= 1; j := 1;
	    end if;
	end while;

	// ensure A[i,j] divides everything in its row and column
	r := 1; s := 1;
	while r le m or s le n do
	    if r le m then
		if A[r,j] mod A[i,j] ne 0 then
		    q := A[r,j] div A[i,j];
		    AddRow(~A,-q,i,r);
		    i := r; r := 1;
		else
		    r +:= 1;
		end if;
	    elif A[i,s] mod A[i,j] ne 0 then
		q := A[i,s] div A[i,j];
		AddColumn(~A,-q,j,s);
		j := s; r := 1; s := 1;
	    else
		s +:= 1;
	    end if;
	end while;

	SwapRows(~A,i,1);
	SwapColumns(~A,j,1);

	// clear out entries below and to the right of A[1,1]
	for i in [2..m] do
	    q := A[i,1] div A[1,1];
	    AddRow(~A,-q,1,i);
	end for;
	for j in [2..n] do
	    q := A[1,j] div A[1,1];
	    AddColumn(~A,-q,1,j);
	end for;

	if A[1,1] lt 0 then
	    A[1,1] := -A[1,1];
	end if;

	C := Submatrix(A, 2, 2, m-1, n-1);
	Diagonalise(~C);
	InsertBlock(~A,C,2,2);
    end if;
end procedure;

procedure SmithNormalForm(~A)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    Diagonalise(~A);
    r := 1;
    while r le Min(m,n) do
	if A[r,r] eq 0 then
	    break;
	else
	    r +:= 1;
	end if;
    end while;
    r -:= 1;

    // enforce divisibility condition
    for i in [1..r] do
	for j in [i+1..r] do
	    if A[j,j] mod A[i,i] ne 0 then
		d := Gcd(A[i,i],A[j,j]);
		l := A[i,i]*A[j,j] div d;
		A[i,i] := d; A[j,j] := l;
	    end if;
	end for;
    end for;
end procedure;
