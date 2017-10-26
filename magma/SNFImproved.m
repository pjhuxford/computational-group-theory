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

	return Sqrt(Min(prod, max^m));
    end if;
end function;

procedure DiagonaliseMod(~A,d)
    Z := IntegerRing(); Zd := IntegerRing(d);
    m := NumberOfRows(A); n := NumberOfColumns(A);

    if not IsZero(A) and m ne 0 and n ne 0 then
	A := ChangeRing(A, Zd);
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
		if (Z ! A[r,j]) mod (Z ! A[i,j]) ne 0 then
		    q := (Z ! A[r,j]) div (Z ! A[i,j]);
		    AddRow(~A,-q,i,r);
		    i := r; r := 1;
		else
		    r +:= 1;
		end if;
	    elif (Z ! A[i,s]) mod (Z ! A[i,j]) ne 0 then
		q := (Z ! A[i,s]) div (Z ! A[i,j]);
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
	    q := (Z ! A[i,1]) div (Z ! A[1,1]);
	    AddRow(~A,-q,1,i);
	end for;
	for j in [2..n] do
	    q := (Z ! A[1,j]) div (Z ! A[1,1]);
	    AddColumn(~A,-q,1,j);
	end for;

	A := ChangeRing(A, Z);

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

	submatrices := {};
	rank := 0;

	for p in primes do
	    // row reduction modulo p
	    // find the p-rank r of A
	    // and an rxr submatrix with full rank (unless we have enough)
	    rows := [1..m];
	    columns := [];

	    Zp := IntegerRing(p);
	    Ap := ChangeRing(A, Zp);

	    i := 1; j := 1;
	    while i le m and j le n do
		k := i;
		while k le m do
		    if Ap[k,j] eq 0 then
			k +:= 1;
		    else
			break;
		    end if;
		end while;

		if k le m then
		    Append(~columns, j);
		    if k ne i then
			SwapRows(~Ap, i, k);
			tmp := rows[k];
			rows[k] := rows[i];
			rows[i] := tmp;
		    end if;
		    c := -1/Ap[i,j];
		    for l in [i+1..m] do
			AddRow(~Ap, Ap[l,j]*c, i, l);
		    end for;
		    i +:= 1;
		end if;
		j +:= 1;
	    end while;

	    r := i - 1;
	    rows := rows[[1..r]];
	    Sort(~rows);
	    if r gt rank then
		rank := r;
		submatrices := {[rows, columns]};
	    elif r eq rank and #submatrices lt dnumber then
		Include(~submatrices, [rows, columns]);
	    end if;
	end for;

	// compute determinants of submatrices
	dets := [];
	for indices in submatrices do
	    rows := indices[1]; columns := indices[2];
	    detsprimes := [];
	    for p in primes do
		Z := IntegerRing(); Zp := IntegerRing(p);
		Ap := ChangeRing(A, Zp);
		Append(~detsprimes, Z ! Determinant(Submatrix(Ap, rows, columns)));
	    end for;
	    sol := CRT(detsprimes, primes);
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
	s := 1;
	while s le Min(m,n) do
	    if A[s,s] eq 0 then
		break;
	    else
		s +:= 1;
	    end if;
	end while;
	s -:= 1;
	// enforce divisibility condition
	for i in [1..s] do
	    for j in [i+1..s] do
		if A[j,j] mod A[i,i] ne 0 then
		    g := Gcd(A[i,i],A[j,j]);
		    l := A[i,i]*A[j,j] div g;
		    A[i,i] := g; A[j,j] := l;
		end if;
	    end for;
	end for;

	// Recover SNF of A
	for i in [s+1..rank] do
	    A[i,i] := d;
	end for;
    end if;
end procedure;
