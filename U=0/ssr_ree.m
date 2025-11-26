function [E_pSSR,E_nSSR] = ssr_ree(RDM)
    % Trim all off-diagonals
    nRDM = diag(diag(RDM));
    % Restore the spin-flip terms
    nRDM(7,10) = RDM(7,10);
    nRDM(10,7) = RDM(10,7);
    % Rotate the spin-flip terms
    aa = nRDM(7,7);
    bb = nRDM(10,10);
    ab = nRDM(7,10);
    ba = nRDM(10,7);
    AB = [aa,ab;ba,bb];
    rot = [1,-1; 1,1]./sqrt(2);
    new = rot*AB*rot';
    nRDM(7,7) = new(1,1);
    nRDM(10,10) = new(2,2);
    pRDM = nRDM;
    % Restore the pair-hopping terms
    pRDM(6,11) = RDM(6,11);
    pRDM(11,6) = RDM(11,6);
    % Rotate the pair-hopping terms
    aa = pRDM(6,6);
    bb = pRDM(11,11);
    ab = pRDM(6,11);
    ba = pRDM(11,6);
    AB = [aa,ab;ba,bb];
    rot = [1,-1; 1,1]./sqrt(2);
    new = rot*AB*rot';
    pRDM(6,6) = new(1,1);
    pRDM(11,11) = new(2,2);
    p = diag(nRDM);
    q = diag(nRDM);
    p(7)  = max(nRDM(7,7),nRDM(10,10));
    p(10) = min(nRDM(7,7),nRDM(10,10));
    a1 = p(7)+p(10)+p(4)+p(13);
    b1 = a1^2-(p(4)-p(13))^2;
    c1 = (p(7)-p(10))*a1;
    d1 = (p(4)+p(13))^2*(p(7)-p(10))^2+8*p(4)*p(13)*(2*p(4)*p(13)+(p(4)+p(13))*(p(7)-p(10))+2*p(7)*p(10));
    q(7) = (b1+c1+sqrt(d1))/(4*(a1-p(10)));
    q(10) = (b1-c1-sqrt(d1))/(4*(a1-p(7)));
    q(4) = p(4) + (p(7)+p(10)-q(7)-q(10))/2;
    q(13) = p(13) + (p(7)+p(10)-q(7)-q(10))/2;
    % P-SSR computation, general formula
    p(6)  = max(pRDM(6,6),pRDM(11,11));
    p(11) = min(pRDM(6,6),pRDM(11,11));
    a2 = p(6)+p(11)+p(1)+p(16);
    b2 = a2^2-(p(1)-p(16))^2;
    c2 = (p(6)-p(11))*a2;
    d2 = (p(1)+p(16))^2*(p(6)-p(11))^2+8*p(1)*p(16)*(2*p(1)*p(16)+(p(1)+p(16))*(p(6)-p(11))+2*p(6)*p(11));
    q(6) = (b2+c2+sqrt(d2))/(4*(a2-p(11)));
    q(11) = (b2-c2-sqrt(d2))/(4*(a2-p(6)));
    q(1) = p(1) + (p(6)+p(11)-q(6)-q(11))/2;
    q(16) = p(16) + (p(6)+p(11)-q(6)-q(11))/2;
    % Global singlet assertion
    if all(abs(RDM(4,4)-RDM(13,13))<1e-25)
        %% N-SSR computation, assuming global singlet
        t = p(7);
        r = p(10) + p(4) + p(13);
        E_nSSR = 0;
        if r < t
            E_nSSR = E_nSSR + r * log2(2*r/(r+t)) + t * log2(2*t/(r+t));
        end
    else
        %% N-SSR computation, general formula
        E_nSSR = 0; indices = [4,7,10,13];
        for i=indices
            if(p(i)==0)
                continue
            end
            E_nSSR = E_nSSR + p(i)*log2(p(i)/q(i));
        end
    end
    % Particle-hole symmetry assertion
    if all(abs(RDM(1,1)-RDM(16,16))<1e-25)
        %% P-SSR computation, assuming particle-hole symmetry
        T = p(6);
        R = p(11) + p(1) + p(16);
        E_pSSR = 0;
        if r < t
            E_pSSR = E_pSSR + R * log2(2*R/(R+T)) + T * log2(2*T/(R+T));
        end
    else
        %% P-SSR computation, general formula
        E_pSSR = 0; indices = [1,6,11,16];
        for i=indices
            if(p(i)==0)
                continue
            end
            E_pSSR = E_pSSR + p(i)*log2(p(i)/q(i));
        end
    end
    % PPT CRITERION in each charge sector
    [N0,~,N2] = sym_negativity(RDM);
    if N0 < 1e-8
        E_nSSR = 0;
    end
    if N2 < 1e-8
        E_pSSR = E_nSSR;
    else
        E_pSSR = E_nSSR + E_pSSR;
    end
end