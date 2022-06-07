function omegahat = ml(signal)

%Type cast it to double.
signal = double(signal);

%--------------------------------------------------------
% Error condition checks
[T,N] = size(signal);

if any([T N]==0)
  error('ml: zero size data not allowed.');
end

if T==1
   signal = signal(:);
   [T,N] = size(signal);
end

vect = []
end