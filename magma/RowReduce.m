RowReduce := procedure(~A)
    m := NumberOfRows(A); n := NumberOfColumns(A);
    i := 1; j := 1;

    while i le m and j le n do
        AbsEntries := [ Abs(A[k,j]) : k in [i..m] ];
	Indices := [i..m];
	ParallelSort(~AbsEntries, ~Indices);
	while true do

	    if AbsEntries[#AbsEntries] eq 0 then
		j +:= 1;
		break;
	    elif i eq m or AbsEntries[#AbsEntries-1] eq 0 then
		k := Indices[#Indices];
		if k ne i then
		    SwapRows(~A,i,k);
		end if;
		if A[i,j] lt 0 then
		    MultiplyRow(~A,-1,i);
		end if;

		for l in [1..i-1] do
		    q := A[l,j] div A[i,j];
		    AddRow(~A,-q,i,l);
		end for;
		i +:= 1; j +:= 1;
		break;

	    else
		k := Indices[#Indices-1]; l := Indices[#Indices];
		q := A[l,j] div A[k,j];
		AddRow(~A,-q,k,l);
		AbsEntries[#AbsEntries] := Abs(A[l,j]);
		ParallelSort(~AbsEntries, ~Indices);
	    end if;
	end while;
    end while;
end procedure;
