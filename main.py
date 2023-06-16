import streamlit as st
import pandas as pd
from sympy import *
import sympy as sy
import scipy as sp
import numpy as np 

x,y,z,t, c = symbols('x y z t c', real=True)
f,g,h = symbols('f g h', cls=Function)
lam = symbols("lambda", real=True)
#init_printing(use_unicode=True)


titre = st.title(" Optimum : MathematiCa")
dessous = st.text("""Cette application vous permettra de résoudre les problèmes\nd'optimisation avec ou sans contraintes""") 

fonction = st.text_input("Entrez la fonction : ")  # Variable pour demander a l'utilisateur de saisir les fonctions
fg = fonction.split(';') #Ici on va convertir fonction en liste pour permettre a l'utilisateur de saisir la fonction et la contrainte
systeme = fonction.split(',')
if ";" in fonction: 
   f,g = fg[0],fg[1] #Ici on sépare chaque fonction de la liste
   f,g = sympify(f, evaluate=False),sympify(g, evaluate=False) #Ici on convertit les deux fonctions en fonction sympy 
   
   #Ici on affiche les deux fonctions
   fon, gon = latex(f), latex(g) # Creer juste pour permettre un affichage en latex
   
   st.subheader("La fonction à optimiser : ")
   st.latex(fr"""\begin{{cases}}f(x,y):{fon}\\ \;s/c \\g(x,y):{gon}=0 \end{{cases}}""")
    
   # Calcul du gradient et de l matrice hessienne
   # Chercher a afficher les matrices 
   # Ici on convertit d'abord les matrix en expression latex
   lagrange = f+(lam*g) # Fonction lagrangienne
   f_lag = latex(lagrange) # Pour afficher la fonction lagrange en latex
   grad_lag = [diff(lagrange,c) for c in ['x','y']]#,'z']]
   kkt_eqs = grad_lag+[g]
   
   #Derivées avec contraintes
   Lx = diff(lagrange, 'x') #Dérivées premiere de L par rapport  à x
   Lx2 = diff(Lx, 'x') #Dérivées seconde de L par rapport à x
   Lxy = diff(Lx, 'y') ##Dérivées seconde croisé de L par rapport à y
   
   
   Ly = diff(lagrange, 'y') #Dérivées premiere de L par rapport  à y
   Ly2 = diff(Ly, 'y') #Dérivées seconde de L par rapport à y
   Lyx = diff(Ly, 'x') ##Dérivées seconde croisé de L par rapport à x
   
   #Dérivées de la fonction g
   gx = diff(g,'x')
   gy = diff(g,'y')
  
 
   hess = Matrix([[Lx2,Lxy,gx],[Lyx, Ly2,gy],[gx,gy,0]])
   hessien = latex(hess)
   

   st.markdown(fr"Ecrivons la fonction lagrangienne : $L(x,y,\lambda) = f(x)+\lambda g(x)$")
   st.markdown(fr"""La fonction de Lagrange  est donc : $L(x,y,\lambda) = {f_lag}$""")
   
   # A.1
   
   #choix = st.selectbox("Condition Premier Ordre", ("oui", "non"),index=1)
   #if choix == "oui":
   st.subheader("Condition de premier ordre : ")
   st.markdown(r"""$\begin{cases} L'_{x}(x,y,\lambda)=0 \\L'_{y}(x,y,\lambda)=0 \\L'_{\lambda}(x,y,\lambda)=0\end{cases}\Rightarrow 
                  \begin{cases}f'_{x}(x,y)+\lambda g'_{x}(x,y)=0 \\f'_{y}(x,y)+\lambda g'_{y}(x,y)=0 \\ g(x,y)=0\end{cases}$""")
   # A.2
   
   st.subheader("""Le systeme à résoudre est le suivant : """)
      
   #      Ici on sait que le systeme à résoudre se trouve dans la variable kkt_eqs
   #      donc on sait que cette variable est une liste de trois éléments qu'on peut facilement
   #      afficher sous forme de systeme d'equation avec latex.
   #      Pour cela on isole éléments par éléménts dans les variable k0,k1 et k2
   #      puis on convertit respectivement ces variable (k0,k1,k2) en latex l0,l1 et l2
   #      pour afficher correctement le systeme à résoudre
      
   k0,k1,k2 = kkt_eqs[0],kkt_eqs[1],kkt_eqs[2]
   l0,l1,l2 = latex(k0), latex(k1), latex(k2)
   st.latex(fr"""S(x,y, \lambda) = \begin{{cases}}{l0} = 0  \\ {l1} = 0 \\ {l2} = 0 \end{{cases}}""")
      
         
   st.write("La matrice hessienne est : ")
   st.latex(fr"""D²(x,y,\lambda) = {hessien}""")   
      
   st.write("Les points stationnaires sont : ")
   point = solve(kkt_eqs, ['x','y',lam], dict=True) #,'z',lam], dict=True)
   solutions = []
   for elt in point:
      elt = latex(elt)
      solutions.append(elt)
      st.latex(fr"""(x,y,\lambda) \in {elt}""")
      
   st.subheader("Condition de second ordre : ") 
   
   determinant = det(hess)
   
   #st.write(determinant)
   #st.write(solutions)
   #print(solutions)
   #print()
   
   #st.code(determinant)
   #v = determinant.subs({'x':2,'y':3,'lam':-108})
   #st.write("# Delta")
   #st.write(v)
   
   det_finale = []
   for sol in point:
      s = determinant.subs(sol)
      det_finale.append(s)
   
   total = sum(det_finale)
   
   
   
   
   
  
   if total < 0:
      #st.write("Le determinant de notre matrice est : ")
      #st.write(total)
      de = latex(determinant)
      st.markdown(fr"""Le determinant de la matrice  $\Delta = {de} < 0$""")
      #st.markdown(fr"""$\Delta = {determinant} < 0$""")
      st.markdown(fr"""Alors le point candidat est un minimum""")
   else:
      #st.write("Le determinant de notre matrice est : ")
      #st.write(total)
      de = latex(determinant)
      st.markdown(fr"""Le determinant de la matrice  $\Delta = {de} > 0$""")
      #st.markdown(fr"""$\Delta = {determinant} > 0$""")
      st.markdown(fr"""Alors le point candidat est un maximum""")
      
elif "," in fonction:      # Nouvelle focntionnalité : Résoudre les système d'équations
   st.subheader("Résolution du système : ")
   equa1, equa2 = systeme[0],systeme[1]
   equation1= sy.sympify(equa1,evaluate=False) # Conversion en fonction sympy
   equation2 = sy.sympify(equa2,evaluate=False) # Conversion en fonction sympy
   latex_e1, latex_e2 = latex(equation1),latex(equation2) # Transformation en latex pour un meilleur affichage
   st.latex(fr"""S : \begin{{cases}}{latex_e1}=0\\ \\{latex_e2}=0 \end{{cases}}""")
   # Résolution du systeme
   eq1 = sy.Eq(equation1, 0) # Transformation en donnée pour résoudre le systeme
   eq2 = sy.Eq(equation2, 0) # Transformation en donnée pour résoudre le systeme
   solution = sy.solve((eq1,eq2),('x','y'))
   st.subheader("Les solutions sont : ")
   for sol in solution:
      s = latex(sol)
      st.latex(fr"""S(x,y) = {s}""")
      
   
   
   
   

else:
   try:
      #fx = simplify(fonction) #ici la fonction sera simplifier 
      f = sympify(fonction) #, evaluate=False) #ici la fonction n'a pas de simplification
   except ValueError:
      st.warning("Veuillez saisir une fonction ") #Test si il n'y a aucune fonction saisie
   else:
      # Optimisation à une variable sans contraintes
      if 'x' in fonction and 'y' not in fonction :
         st.write("# La fonction à optimiser : ")
         st.write("Optimisation sans contrainte")
         
         fonc = latex(f) #Ici on transforme notre fonction en code latex pour mieux l'afficher dans le programme
         st.latex(fr"""f(x) = {fonc}""")
         #st.write(f)
         
         st.subheader("Condition de premier ordre : ")
         st.markdown(r"""L'equation $f'(x)=0$ doit admettre au moins une solution""") #Ici on affiche seulement f'(x) = 0

         #Pour mieux afficher les derivees en latex
         fx = diff(f,'x')
         fx2 = diff(fx, 'x')
         dx = latex(fx)  # Transformation en latex pour mieux afficher la derivee dans le programme
      
         st.markdown(fr"""La dérivée de notre fonction est : """) #{fx} \Leftrightarrow f'(x) = 0$""")            
         st.latex(fr"""f'_{'x'} = {dx}""")      
         st.markdown(fr"""Résolution de l'équation $f'(x) = 0$""")
         equation = solve(fx, 'x')
      
         st.text("Solution(s) de l'équation : ")
         i = 0
         S = [] #Liste qui regroupera toutes les solutions
         for solution in equation:
            solution = latex(solution)
            st.latex(fr"""x_{i} = {solution}""")
            S.append(solution)
            i += 1
         
         st.subheader("Condition de second ordre : ")
         fx2 = diff(fx, 'x') #Derivee seconde de f 
            #<> 
         st.markdown(r"""La dérivée seconde de notre fonction est : """)
         dx2 = latex(fx2)
            
         st.latex(fr"""f''_{'x'} = {dx2}""")
         for xi in S:
            fxi = fx2.subs('x',xi)
            if fxi < 0 :
               st.markdown(fr"""$f''(x_{"0"}) \Leftrightarrow f''({xi}) = {fxi} < 0 \\$ Donc le point $x_{"0"} = {xi}$ est un maximum et la fonction est concave""")      
            else:
               st.markdown(fr"""$f''(x_{"0"}) \Leftrightarrow f''({xi}) = {fxi} >  0 \\$ Donc le point $x_{"0"} = {xi}$ est un minimum et la fonction est convexe""")
            i += 1   
      else:
         # optimisation à deux variables sans contraintes
         st.write("# La fonction à optimiser : ")
         st.write("Optimisation sans contrainte")
         fonc = latex(f)
         st.latex(fr"""f(x,y) = {fonc}""")
          
         st.subheader("Condition du premier ordre ")
         
         #Calcul des derivees de notre fonction
         fx = diff(f, 'x') # Derivee premiere par rapport a x
         fx2 = diff(fx, 'x') # Derivee seconde par rapport a x
         fxy = diff(fx, 'y') # Derivees croisee de fx par rapport a y
         
         fy = diff(f, 'y') # Derivee premiere par rapport a y
         fy2 = diff(fy, 'y') # Derivee seconde par rapport a y
         fyx = diff(fy, 'x') # Derivees croisee de fy par rapport a x
         
         #Affichage en latex pour les derivees
         dx = latex(fx)
         dx2 = latex(fx2)
         dxy = latex(fxy)
         #  ici on met ça pour mieux les afficher en latex dans le programme
         
         dy = latex(fy)
         dy2 = latex(fy2)
         dyx = latex(fyx)
         hess = Matrix([[fx2,fxy],[fyx,fy2]])
         hessien = latex(hess)
         
         # Ici on met ça pour afficher le gradient
         grad = Matrix([[fx],[fy]]) 
         gradient = latex(grad) # Ici c'est pour afficher le gradient en latex
         
         # Ici on affiche la condition du premier ordre en latex  
         st.text("Les dérivées de notre fonction sont : ")
         st.latex(fr"""f'_{{x}} = {dx} \quad f''_{{x}} = {dx2} \quad f''_{{xy}} = {dxy} \\
                  f'_{{y}} = {dy} \quad f''_{{y}} = {dy2} \quad f''_{{yx}} = {dyx}""")
         
         
         st.text("Pour que la fonction admet un extremum :\nIl faut que le systeme ci-dessous admet au moins une solution")          
         st.latex(fr"""\begin{{cases}} f'_{x}(x,y)=0 \\f'_{y}(x,y)=0 \end{{cases}}\Rightarrow 
               \begin{{cases}}{dx} = 0\\ {dy} = 0\end{{cases}}""")
         
         st.subheader("Le gradient de notre fonction est : ")
         st.latex(fr"""D(x,y) = {gradient}""")
         
         st.subheader("La matrice hessienne est : ")
         st.latex(fr"""H²(x,y) = {hessien}""")
         
         # resolution du systeme d'equation de la condition du second ordre
         eq1 = Eq(fx,0)
         eq2 = Eq(fy,0)
         #solution = solve(grad,[x,y], dict=True)
         
         solution = solve((fx,fy), ("x","y"), dict=True)
         print(solution)
         
         solutions = solution[0]# Ici solution est une liste , on fait solutions[0] pour qu'elle soit un dictionnaire
        
         
        
        
         st.subheader("Les solutions du systeme sont : ")
         #st.latex(fr"""S(x,y) = {solution}""")
         for sol in solution:
            s = latex(sol)
            st.latex(fr"""S(x,y) = {s}""")
            
       
         
         #   st.latex(fr"""(x,y) \in {sol}""")
         #   point_solution.append(sol)
         #st.write(point_solution)
         
         s = [] # Liste qui contiendra les solutions
         for valeur in solutions.values():
            s.append(valeur)
         x0,y0 = s[0],s[1]
         x0,y0 = int(x0),int(y0)
         
            
         
      
         st.subheader("Condition du second ordre ")
         
         # METHODES 
         r0 = fx2.subs({'x':x0,'y':y0})
         s0 = fxy.subs({'x':x0,'y':y0})
         t0 = fy2.subs({'x':x0,'y':y0})
         
         s = r0*t0-(s0*s0)
         s2 = r0+t0
         st.markdown(r"""La nature de notre point depend du signe de $\Delta f$""")
         if s > 0 and s2 > 0: 
            st.markdown(fr"""Comme $\Delta f$ = {s} est strictement superieur à 0 et que $r+t = {s2} > 0$, Alors le point ({x0};{y0}) est un minimum local""")
         elif s > 0 and s2 < 0:
             st.markdown(fr"""Comme $\Delta f$ = {s} est strictement superieur à 0 et que $r+t = {s2} < 0$, Alors le point ({x0};{y0}) est un maximum local""")
         elif s < 0 :
            st.markdown(fr"""Comme $\Delta f$ = {s} est strictement inferieur à 0, Alors le point ({x0};{y0}) est un point selle""")
         else:
            st.markdown(fr"""Comme $\Delta f$ = {s} est egale à 0, On ne peut conclure""")
        
            
         
            
         
         
bas = st.caption("Développez par Dannys")