procedure Diagonalise(~A, ~max, ~overtime, maxtime)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    if not IsZero(A) and m ne 0 and n ne 0 then
	// find min nonzero modulus entry A[i,j]
	// store max modulus entry
	r := 1; s := 1;
	while A[r,s] eq 0 do
	    if s lt n then
		s +:= 1;
	    else
		r +:= 1; s := 1;
	    end if;
	end while;
	i := r; j := s;
	max := Abs(A[r,s]);
	min := max;
	while r le n do
	    if A[r,s] ne 0 and Abs(A[r,s]) lt min then
		min := Abs(A[r,s]);
		i := r; j := s;
	    elif Abs(A[r,s]) gt max then
		max := Abs(A[r,s]);
	    end if;
	    if s lt n then
		s +:= 1;
	    else
		r +:= 1; s := 1;
	    end if;
	end while;

	// ensure A[i,j] divides everything in its row and column
	r := 1; s := 1;
	while r le m or s le n do
	    if Cputime() gt maxtime then
		overtime := true;
		print "overtime";
		break;
	    else
		if r le m then
		    if A[r,j] mod A[i,j] ne 0 then
			q := A[r,j] div A[i,j];
			AddRow(~A,-q,i,r);
			i := r; r := 1;
			for k in [1..n] do
			    if Abs(A[i,k]) gt max then
				max := Abs(A[i,k]);
			    end if;
			end for;
		    else
			r +:= 1;
		    end if;
		elif A[i,s] mod A[i,j] ne 0 then
		    q := A[i,s] div A[i,j];
		    AddColumn(~A,-q,j,s);
		    j := s; r := 1; s := 1;
		    for k in [1..m] do
			if Abs(A[k,j]) gt max then
			    max := Abs(A[k,j]);
			end if;
		    end for;
		else
		    s +:= 1;
		end if;
	    end if;
	end while;

	if not overtime then
	    SwapRows(~A,i,1);
	    SwapColumns(~A,j,1);

	    // clear out entries below and to the right of A[1,1]
	    for i in [2..m] do
		q := A[i,1] div A[1,1];
		AddRow(~A,-q,1,i);
		for k in [1..n] do
		    if Abs(A[i,k]) gt max then
			max := Abs(A[i,k]);
		    end if;
		end for;
	    end for;
	    for j in [2..n] do
		q := A[1,j] div A[1,1];
		AddColumn(~A,-q,1,j);
		for k in [1..m] do
		    if Abs(A[k,j]) gt max then
			max := Abs(A[k,j]);
		    end if;
		end for;
	    end for;

	    if A[1,1] lt 0 then
		A[1,1] := -A[1,1];
	    end if;

	    C := Submatrix(A, 2, 2, m-1, n-1);
	    submax := 0;
	    Diagonalise(~C, ~submax, ~overtime, maxtime);
	    InsertBlock(~A,C,2,2);
	    if submax gt max then
		max := submax;
	    end if;
	end if;
    end if;
end procedure;

procedure SmithNormalForm(~A,~max : maxtime := 1800)
    overtime := false;
    m := NumberOfRows(A); n := NumberOfColumns(A);
    Diagonalise(~A,~max,~overtime,maxtime);
    if not overtime then
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
    end if;
end procedure;
