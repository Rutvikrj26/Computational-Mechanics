function [out] = force_check(node_no)
  z = load('nodecoordinates_square_400.txt');
  h = 0.2;
  out = 0;
  if z(node_no,3)>=(-0.5+h) && z(node_no,3)<=(0.5-h)
    out = 1;
  endif
endfunction