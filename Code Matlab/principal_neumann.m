% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
Fichiers=["geomCarreh04.msh" "geomCarreh02.msh" "geomCarreh01.msh" "geomCarreh005.msh"];
h=[1/0.4 1/0.2 1/0.1 1/0.05];
erreurL2=zeros(1,4);
erreurH1=zeros(1,4);

for k=1:4
    nom_maillage = Fichiers(k)
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

    % ----------------------
    % calcul des matrices EF
    % ----------------------

    % declarations
    % ------------
    KK = sparse(Nbpt,Nbpt); % matrice de rigidite
    MM = sparse(Nbpt,Nbpt); % matrice de rigidite
    LL = zeros(Nbpt,1);     % vecteur second membre

    % boucle sur les triangles
    % ------------------------
    for l=1:Nbtri
      % Coordonnees des sommets du triangles
      % A COMPLETER
      S1=Coorneu(Numtri(l,1),:);
      S2=Coorneu(Numtri(l,2),:);
      S3=Coorneu(Numtri(l,3),:);
      % calcul des matrices elementaires du triangle l 

       Kel=matK_elem(S1, S2, S3);
       Mel=matM_elem(S1, S2, S3);

      % On fait l'assemmblage de la matrice globale et du second membre
      % A COMPLETER
      for i=1:3
          for j=1:3
            MM(Numtri(l,i),Numtri(l,j))=MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
            KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
          end
      end

    end % for l

    % Calcul du second membre L
    % -------------------------
    % A COMPLETER
    % utiliser la routine f.m
    FF = zeros(Nbpt,1);
    for i= 1:Nbpt
        FF(i) = f(Coorneu(i,1), Coorneu(i,2));
    end;    
    LL = MM*FF;

    % inversion
    % ----------
    UU = (MM+KK)\LL;

    % visualisation
    % -------------
    affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

    validation = 'oui';
    % validation
    % ----------
    if strcmp(validation,'oui')
        UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
        %affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
        % Calcul de l erreur L2
        % A COMPLETER
        L2 = sqrt(transpose(UU_exact - UU)*MM*(UU_exact - UU));
        % Calcul de l erreur H1
        % A COMPLETER
        % attention de bien changer le terme source (dans FF)
        H1 = sqrt(transpose(UU_exact - UU)*KK*(UU_exact - UU));
    end
    normL2_U = sqrt(transpose(UU_exact)*MM*(UU_exact));
    normH1_U = sqrt(transpose(UU_exact)*KK*(UU_exact));
    
    erreurL2(k)=L2/normL2_U;
    erreurH1(k)=H1/normH1_U;
end
%regressionL2brut=polyfit(log(h),log(erreurL2),1);
%regressionH1brut=polyfit(log(h),log(erreurH1),1);

%hregression = linspace(log(h(1)),log(h(4)));
%regressionL2=polyval(regressionL2brut,hregression);
%regressionH1=polyval(regressionH1brut,hregression);

%grid on
%hold on
%plot(log(h),log(erreurL2),'o-');
%plot(log(h),log(erreurH1),'*-');
%plot(hregression,regressionL2,'r--');
%plot(hregression,regressionH1,'b--');
%legend('L2 norm','H1 seminorm',sprintf('L2 norm lin m=%f',regressionL2brut(1)),sprintf('H1 norm lin m=%f',regressionH1brut(1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023

