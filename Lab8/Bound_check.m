function [der_x,der_y] = Bound_check(vect)
  der_x = 0;
  der_y = 0;
  if vect(3)==0
    der_x = 1;
  end
  if vect(4)==0
    der_y = 1;
  end
endfunction