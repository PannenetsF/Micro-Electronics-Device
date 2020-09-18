function y_data = fakelog(y)
  eps = 1e-16;
  y_t = abs(y) + eps; 
  y_data = log10(y_t) .* ((y>0)*2-1);
endfunction
