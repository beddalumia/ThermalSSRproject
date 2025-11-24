function print_basis()
%% Pretty print of two-site basis states, in EDIpack ordering
for state = 0:1:15
   [~,label] = build_ket(state,2);
   fprintf('%d\t',state+1)
   disp(label)
end
end

% FROM ED_SETUP:
% |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
% |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
% lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
