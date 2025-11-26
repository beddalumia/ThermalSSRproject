function [FT,PT,new_indices] = partial_transpose(RDM_ij)
    % No Fermi signs, no phases! Just move around elements for the "bosonic" PT
    [Nrdm,Mrdm] = size(RDM_ij); assert(Nrdm==Mrdm);
    ROT = zeros(Nrdm,Mrdm); PT = ROT;
    % Rotate the RDM_ij to the proper block form
    [ROT,new_indices,signs] = single_site_sectors(RDM_ij);
    % Fermionic tranpose: we need to take account of signs and phases
    Norb = 1;
    Nlso = 4; % 2 sites x 2 spin x 1 orbital
    for i = 1:16
        ket = fliplr(dec2bin(new_indices(i)-1,Nlso));
        k_A = ket(1:2*Norb);
        k_B = ket(2*Norb+1:end);
        P_Ai = count(k_A,'1');
        P_Bi = count(k_B,'1');
        for j = 1:16
            bra = fliplr(dec2bin(new_indices(j)-1,Nlso));
            b_A = bra(1:2*Norb);
            b_B = bra(2*Norb+1:end);
            P_Aj = count(b_A,'1');
            P_Bj = count(b_B,'1');
            % Fermi phase (twisted)
            FROT(new_indices(i),new_indices(j)) = ROT(new_indices(i),new_indices(j)) * (-1)^(0.5*mod(P_Bi+P_Bj,2)+(P_Ai+P_Aj)*(P_Bi+P_Bj) + P_Bi);
        end
    end
    % Finally, the partial transpose on each "left site" block
    N_1site = 4;
    for i = 1:Nrdm/N_1site
       for j = 1:Nrdm/N_1site
          block = ROT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site);
          PT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site) = block';
          block = FROT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site);
          FT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site) = block';
       end
    end
    % Return to EDIpack basis: |iup jup idw jdw> (we want this function to be idempotent!)
    Ftemp = FT; Btemp = PT;
    for i = 1:16
        for j = 1:16
            PT(i,j) = Btemp(new_indices(i),new_indices(j));
            if xor(any(i==signs),any(j==signs))
                % Uncanceled Fermi sign!
                FT(i,j) = -Ftemp(new_indices(i),new_indices(j));
            else
                FT(i,j) = Ftemp(new_indices(i),new_indices(j));
            end
        end
    end
end
function [ROT,new_indices,fermi_signs] = single_site_sectors(RDM_ij)
% Rotate the RDM_ij
% from the |i_up j_up〉⊗ i_dw j_dw〉
% to the   |i_up i_dw〉⊗ j_up j_dw〉
% basis (useful for negativity, partial trace, etc)
% ACHTUNG: take care of all the fermionic signs!
% 1       | •  • 〉⊗ | •  • 〉->    | •  • 〉⊗ | •  • 〉
% 2       | ↑  • 〉⊗ | •  • 〉->    | ↑  • 〉⊗ | •  • 〉
% 3       | •  ↑ 〉⊗ | •  • 〉->    | •  • 〉⊗ | ↑  • 〉
% 4       | ↑  ↑ 〉⊗ | •  • 〉->    | ↑  • 〉⊗ | ↑  • 〉
% 5       | •  • 〉⊗ | ↓  • 〉->    | •  • 〉⊗ | •  ↓ 〉
% 6       | ↑  • 〉⊗ | ↓  • 〉->    | ↑  ↓ 〉⊗ | •  • 〉
% 7       | •  ↑ 〉⊗ | ↓  • 〉-> -1 | •  ↓ 〉⊗ | ↑  • 〉
% 8       | ↑  ↑ 〉⊗ | ↓  • 〉-> -1 | ↑  ↓ 〉⊗ | ↑  • 〉
% 9       | •  • 〉⊗ | •  ↓ 〉->    | •  • 〉⊗ | •  ↓ 〉
% 10      | ↑  • 〉⊗ | •  ↓ 〉->    | ↑  • 〉⊗ | •  ↓ 〉
% 11      | •  ↑ 〉⊗ | •  ↓ 〉->    | •  • 〉⊗ | ↑  ↓ 〉
% 12      | ↑  ↑ 〉⊗ | •  ↓ 〉->    | ↑  • 〉⊗ | ↑  ↓ 〉
% 13      | •  • 〉⊗ | ↓  ↓ 〉->    | •  ↓ 〉⊗ | •  ↓ 〉
% 14      | ↑  • 〉⊗ | ↓  ↓ 〉->    | ↑  ↓ 〉⊗ | •  ↓ 〉
% 15      | •  ↑ 〉⊗ | ↓  ↓ 〉-> -1 | •  ↓ 〉⊗ | ↑  ↓ 〉
% 16      | ↑  ↑ 〉⊗ | ↓  ↓ 〉-> -1 | ↑  ↓ 〉⊗ | ↑  ↓ 〉
%
    fermi_signs = [7,8,15,16];
    [Nrdm,Mrdm] = size(RDM_ij); assert(Nrdm==Mrdm);
    Nlat = 2;
    Norb = 1;
    Nlso = 4; % 2 sites x 2 spin x 1 orbital
    ROT = zeros(Nrdm,Mrdm); %PT = ROT;
    % Rotate the RDM_ij to the proper block form
    for i = 1:Nrdm
       ket = fliplr(dec2bin(i-1,Nlso));
       kup = ket(1:Nlat*Norb);
       kdw = ket(Nlat*Norb+1:end);
       for ik = 1:Nlso
          if(mod(ik,2)==1)
                ket(ik) = kup(round((ik+1)/2));
          else
                ket(ik) = kdw(round(ik/2));
          end
       end
       newI = bin2dec(fliplr(ket))+1;
       new_indices(i) = newI;
       for j = 1:Nrdm
             bra = fliplr(dec2bin(j-1,Nlso));
             bup = bra(1:Nlat*Norb);
             bdw = bra(Nlat*Norb+1:end);
             for ik = 1:Nlso
                if(mod(ik,2)==1)
                   bra(ik) = bup(round((ik+1)/2));
                else
                   bra(ik) = bdw(round(ik/2));
                end
             end
             newJ = bin2dec(fliplr(bra))+1;
             % (i,j) ---> (newI,newJ)
             if xor(any(i==fermi_signs),any(j==fermi_signs))
                % Uncanceled Fermi sign!
                ROT(newI,newJ) = -RDM_ij(i,j);
            else
                ROT(newI,newJ) = RDM_ij(i,j);
            end
       end
    end
end
