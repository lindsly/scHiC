function out = mylog2_neg_inf(in)
%   find(isnan(in))
  out = log2(in+.5)+1;
  out(~in) = 0;
end
