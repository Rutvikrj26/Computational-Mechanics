function bound_yesno = Bound_yesno(vector)
  bound_yesno = 0;
  if nnz(vector) == 4
    bound_yesno = 1;
  end
endfunction