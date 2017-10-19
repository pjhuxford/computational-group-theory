Diagonalise := procedure(~A)
    if not IsZero(A) then
	m := NumberOfRows(A); n := NumberOfColumns(A);
	i := 1; j := 1;
	while A[i,j] eq 0 do
	    if j lt n then
		j +:= 1;
	    else
		i +:= 1;
		j := 1;
	    end if;
	end while;

	r := 1; s := 1;
	while r le m and s le n do
	    if r le m then
		if A[r,j] mod A[i,j] ne 0 then
		    q := A[r,j] div A[i,j];
		    AddRow(~A,-q,i,r);
		    i := r;
		    r := 1;
		else
		    r +:= 1;
		end if;
	    elif A[i,s] mod A[i,j] ne 0 then
		q := A[i,s] div A[i,j];
		AddColumn(~A,-q,j,s);
		j := s;
		r := 1;
		s := 1;
	    else
		s +:= 1;
	    end if;
	end while;
    end if;
end procedure;
