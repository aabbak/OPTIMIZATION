function bin = binaryTransfer(x)
% BINARYTRANSFER  Map a continuous search vector to a 0/1 selection vector.
% Used by every optimiser so the continuous metaheuristics operate on the same
% binary reference-point-selection problem.
bin = mod(round(mod(abs(x),2)),2);
end
