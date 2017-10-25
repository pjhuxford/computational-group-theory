function hadamard(A)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    if m gt n then
	return hadamard(Transpose(A));
    else
	prod := 1;
	for i in [1..m] do
	    sum := 0;
	    for j in [1..n] do
		sum +:= A[i,j]^2;
	    end for;
	    if sum ne 0 then
		prod *:= sum;
	    end if;
	end for;

	max := 0;
	for j in [1..n] do
	    sum := 0;
	    for i in [1..m] do
		sum +:= A[i,j]^2;
	    end for;
	    max := Max(max,sum);
	end for;

	return Min(Sqrt(prod), Sqrt(max^m));
    end if;
end function;

procedure DiagonaliseMod(~A,d)
    Z := IntegerRing();
    Z_d := IntegerRing(d);
    A := ChangeRing(ChangeRing(A,Z_d),Z);
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
		    A := ChangeRing(A,Z_d);
		    AddRow(~A,-q,i,r);
		    A := ChangeRing(A,Z);
		    i := r; r := 1;
		else
		    r +:= 1;
		end if;
	    elif A[i,s] mod A[i,j] ne 0 then
		q := A[i,s] div A[i,j];
		A := ChangeRing(A,Z_d);
		AddColumn(~A,-q,j,s);
		A := ChangeRing(A,Z);
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
	    A := ChangeRing(A,Z_d);
	    AddRow(~A,-q,1,i);
	    A := ChangeRing(A,Z);
	end for;
	for j in [2..n] do
	    q := A[1,j] div A[1,1];
	    A := ChangeRing(A,Z_d);
	    AddColumn(~A,-q,1,j);
	    A := ChangeRing(A,Z);
	end for;

	if A[1,1] lt 0 then
	    A[1,1] := -A[1,1];
	end if;

	C := Submatrix(A, 2, 2, m-1, n-1);
	DiagonaliseMod(~C,d);
	InsertBlock(~A,C,2,2);
    end if;
end procedure;

// Computes the SNF of A using Modular techniques
// primestart determines where to begin looking for primes to do calculations
// dnumber determines the number of nonzero rxr determinants we compute the gcd of,
// where r is the rank of A
procedure SmithNormalFormImproved(~A : primestart := 1000, dnumber := 3)
    X := SmithForm(A);
    m := NumberOfRows(A); n := NumberOfColumns(A);
    if not IsZero(A) and m ne 0 and n ne 0 then
	b := hadamard(A);

	// find sequence of distinct primes whose product exceeds 2*b
	primes := [NextPrime(primestart)];
	prod := primes[1];
	while prod le 2*b do
	    Append(~primes, NextPrime(primes[#primes]));
	    prod *:= primes[#primes];
	end while;

	Submatrices := {};
	rank := 0;

	for p in primes do
	    // find the p-rank r of A
	    // and an rxr submatrix with nonzero determinant (unless we have enough)
	    RowIndices := [1..m];
	    ColumnIndices := [];

	    Z_p := IntegerRing(p);
	    A_p := ChangeRing(A, Z_p);

	    i := 1; j := 1;
	    while i le m and j le n do
		k := i;
		while k le m do
		    if A_p[k,j] eq 0 then
			k +:= 1;
		    else
			break;
		    end if;
		end while;

		if k le m then
		    Append(~ColumnIndices, j);
		    if k ne i then
			SwapRows(~A_p, i, k);
			tmp := RowIndices[k];
			RowIndices[k] := RowIndices[i];
			RowIndices[i] := tmp;
		    end if;
		    c := -1/A_p[i,j];
		    for l in [i+1..m] do
			AddRow(~A_p, A_p[l,j]*c, i, l);
		    end for;
		    i +:= 1;
		end if;
		j +:= 1;
	    end while;

	    r := i - 1;
	    RowIndices := RowIndices[[1..r]];
	    Sort(~RowIndices);
	    if r gt rank then
		rank := r;
		Submatrices := {[RowIndices, ColumnIndices]};
	    elif r eq rank and #Submatrices lt dnumber then
		Include(~Submatrices, [RowIndices, ColumnIndices]);
	    end if;
	end for;

	// compute determinants of Submatrices
	dets := [];
	for indices in Submatrices do
	    dets_p := [];
	    for p in primes do
		Z_p := IntegerRing(p);
		A_p := ChangeRing(A, Z_p);
		Append(~dets_p, IntegerRing() ! Determinant(Submatrix(A, indices[1], indices[2])));
	    end for;
	    sol := CRT(dets_p, primes);
	    altsol := sol - prod;
	    if Abs(altsol) lt sol then
		Append(~dets, altsol);
	    else
		Append(~dets, sol);
	    end if;
	end for;

	d := GCD(dets);

	// Compute the SNF of A modulo d
	DiagonaliseMod(~A,d);
	t := 1;
	while t le Min(m,n) do
	    if A[t,t] eq 0 then
		break;
	    else
		t +:= 1;
	    end if;
	end while;
	t -:= 1;
	// enforce divisibility condition
	for i in [1..t] do
	    for j in [i+1..t] do
		if A[j,j] mod A[i,i] ne 0 then
		    g := Gcd(A[i,i],A[j,j]);
		    l := A[i,i]*A[j,j] div g;
		    A[i,i] := g; A[j,j] := l;
		end if;
	    end for;
	end for;

	// Recover SNF of A
	for i in [t+1..rank] do
	    A[i,i] := d;
	end for;
    end if;
end procedure;
