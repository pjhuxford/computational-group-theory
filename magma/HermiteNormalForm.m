RowReduce := procedure(~A : verbose := false)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    i := 1; j := 1;

    if verbose then
	print "A = ", A;
	print "A is a", m, "x", n, "matrix.";
    end if;

    while i le m and j le n do
	if verbose then
	    print "Looking at row", i, "and column", j;
	end if;
	AbsEntries := [ Abs(A[k,j]) : k in [i..m] ];
	Indices := [i..m];
	ParallelSort(~AbsEntries, ~Indices);
	while true do
	    if verbose then
		print AbsEntries, Indices;
	    end if;

	    if AbsEntries[#AbsEntries] eq 0 then
		if verbose then
		    print "only zeros down here";
		end if;
		j +:= 1;
		break;
	    elif i eq m or AbsEntries[#AbsEntries-1] eq 0 then
		k := Indices[#Indices];
		if verbose then
		    print "unique nonzero thingy";
		    print "k =", k;
		end if;
		if k ne i then
		    SwapRows(~A,i,k);
		    if verbose then
			print "Swapping rows", i, k;
			print A;
		    end if;
		end if;
		if A[i,j] lt 0 then
		    MultiplyRow(~A,-1,i);
		    if verbose then
			print "Negating row", i;
			print A;
		    end if;
		end if;

		for l in [1..i-1] do
		    q := A[l,j] div A[i,j];
		    AddRow(~A,-q,i,l);
		    if verbose then
			print "Subtracting", q, "of row", i, "from row", l;
			print A;
		    end if;
		end for;
		i +:= 1; j +:= 1;
		if verbose then
		    print "moving on...";
		end if;
		break;

	    else
		k := Indices[#Indices-1]; l := Indices[#Indices];
		if verbose then
		    print "two distinct nonzero entries";
		    print "k =", k, "l =", l;
		end if;
		q := A[l,j] div A[k,j];
		AddRow(~A,-q,k,l);
		if verbose then
		    print "aha! we will subtract", q, "of row", k, "away from row", l;
		    print A;
		end if;
		AbsEntries[#AbsEntries] := Abs(A[l,j]);
		ParallelSort(~AbsEntries, ~Indices);
	    end if;
	end while;
    end while;
    if verbose then
	print "all done";
    end if;
end procedure;

HermiteNormalForm := function(A : verbose := false)
    B := A;
    RowReduce(~B : verbose := verbose);
    return B;
end function;
