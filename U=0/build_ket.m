function [vec,ket] = build_ket(state,Nlat)
  % BUILD_KET : Puts together a bit representation of a Slater determinant
  %             (and, if asked, a pretty string to print it to screen...)
  %  
  %  >> [vec,ket] = build_ket(state)
  %
  %  state :: integer representation of a basis state (bits are occupation numbers)
  %
  %  This depends entirely on the Fock basis conventions choosen in all
  %  ED-based solvers from QcmPlab (and many other codes, actually)
  %
  %  Copyright 2022 Gabriele Bellomia
  %
  Norb = 1;
  for ilat = 1:Nlat
      for ispin = 1:2
          shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
          index = shift+(1:Norb);
          vec(index) = bitget(state,index);
      end
  end
  if nargout>1
      kup = num2str(vec(1:Norb*Nlat));
      kdw = num2str(vec(Norb*Nlat+1:end));
      ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
          '| ',strrep(kdw,'1','↓'), ' 〉'];
      ket = strrep(ket,'0','•');
  end
end

% FROM ED_SETUP:
% |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
% |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
% lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat