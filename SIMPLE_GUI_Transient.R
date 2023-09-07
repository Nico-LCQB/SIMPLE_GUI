cat("\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("%                                                                                       %\n")
cat("%       SIMPLE : SImulation de Mutations doubles dans des Populations de LEvures        %\n")
cat("%                                                                                       %\n")
cat("%                     nicolas Agier - V7.12 (Lang model-19-06-2023)                     %\n")
cat("%                                                                                       %\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("\n")
cat("\n")

# dans ce modele, seule l'une des deux filles issues de chaque division peut porter une mutation.
# Lang model
# dans cette version, la population hypermutatrice transitoire est simulée sur plusieurs générations consécutives
# Comme inclus dans la simulation de deux mutations independantes, la simulation de 2 simultanées est inactivee
# Cette version donne le nombre de mutations ayant eu lieu dans chaque realisation

###########################
# Fonction de simulation  #
###########################

######################################
# fonction 2 mutations independantes #
######################################

# memory.limit(size = 10000)

SIMUI <- function(x,y,z,c,h,m,tgen, ngen) # x=r1, y=r2, z=g, c = cycle, h = size transient pop, m = increase of rates, tgen = generation of transiant hyperm, ngen=nombre de generation

{
    #--------------------------------------
    # création des objets pour le calcul  |
    #--------------------------------------

    # population de départ
    motherA=c(0)
    motherB=c(0)
    motherAB=c(0)
    tgen=sample(1:27, size =  1)
    while (h >= 2^(tgen-1))
    {
        tgen=sample(1:27, size =  1)
        #print(tgen)
    }
    print(paste("les transient mutator sont apparus a la generation :  ",tgen))

    #----------------------------------------
    # boucle de croissance mutation A et B  |
    #----------------------------------------

    # initialiser compteur pour compter le nombre de mutants

        nmutA=0
        nmutB=0
        nmutAB=0

    # Boucle de simulation

        for (i in 1:z)
        {
            # print(paste("generation : ", i, "res : ",sum(motherA)," ",sum(motherB)," ",length(motherAB[motherAB>1])," ",nmutA," ",nmutB," "nmutAB))

############ generation entree de regime hypermutateur    
            if(i==tgen) # dissocier les deux populations  
            {
                # pour le calcul du nombre de mutants à cette generation
                nmutAi = 2*sum(motherA)
                nmutBi = 2*sum(motherB)
                nmutABi = 2*length(motherAB[motherAB>1])

                # separer les populations au niveau des mother # il doit manquer un +1 ou -1 quelque part :-)
                popsizen=length(motherA)-h
                motherAnm = motherA[1:popsizen]
                motherBnm = motherB[1:popsizen]
                motherAhm = motherA[as.numeric(popsizen+1):length(motherA)]
                motherBhm = motherB[as.numeric(popsizen+1):length(motherB)]

                # separer les populations hyper mutatrices et les autres et tirer les mutants
                daughterAnm = sample(0:1,length(motherA)-h,prob = c(1-x,x), replace = T) 
                daughterBnm = sample(0:1,length(motherB)-h,prob = c(1-y,y), replace = T) 
                daughterAhm = sample(0:1,h,prob = c(1-x*m,x*m), replace = T)
                daughterBhm = sample(0:1,h,prob = c(1-y*m,y*m), replace = T)

                # trouver les doubles apparus dans les filles (chez les normaux et les hypermutateurs
                tempABnm=daughterAnm+daughterBnm
                tempABnm=replace(tempABnm, tempABnm>1,10)
                tempABnm=tempABnm+motherAnm+motherBnm # ajout des caractere maternel pour ne pas mettre des simultane dans les sequentiels
                candidaten=tempABnm[tempABnm==10]
                if(length(candidaten)>0)
                  { print(paste(length(candidaten)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }

                tempABhm=daughterAhm+daughterBhm
                tempABhm=replace(tempABhm, tempABhm>1,10)
                tempABhm=tempABhm+motherAhm+motherBhm # ajout pour ne pas mettre des simultane dans les sequentiels
                candidateh=tempABhm[tempABhm==10]
                if(length(candidateh)>0)
                  { print(paste(length(candidateh)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }

                # nettoyage des fichiers temp
                rm(tempABnm)
                rm(tempABhm)
                
                # transfert mother state to daughter state
                daughterAnm = motherAnm + daughterAnm
                daughterBnm = motherBnm + daughterBnm
                daughterAhm = motherAhm + daughterAhm
                daughterBhm = motherBhm + daughterBhm

                # junction of 2 "daughter" lists
                motherAnm = c(motherAnm, daughterAnm)
                motherBnm = c(motherBnm, daughterBnm) 
                motherAhm = c(motherAhm, daughterAhm)
                motherBhm = c(motherBhm, daughterBhm) 

                # correct for mutation in mother already mutated
                motherAnm=replace(motherAnm, motherAnm>0,1)
                motherBnm=replace(motherBnm, motherBnm>0,1)
                motherAhm=replace(motherAhm, motherAhm>0,1)
                motherBhm=replace(motherBhm, motherBhm>0,1)

                # table for dble mutant
                motherABnm=motherAnm+motherBnm
                motherABhm=motherAhm+motherBhm

                # nettoyage des fichiers temp
                rm(daughterAnm)        
                rm(daughterBnm)
                rm(daughterAhm)        
                rm(daughterBhm)

                # calcul du nombre de mutant à cette generation # ICI ON ADDITIONNE LES DEUX OU REMERGE LES DEUX POPULATIONS
                nmutAf = sum(motherAnm)+sum(motherAhm)
                nmutBf = sum(motherBnm)+sum(motherBhm)
                nmutABf = length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
                nmutA=nmutA+nmutAf-nmutAi
                nmutB=nmutB+nmutBf-nmutBi
                nmutAB=nmutAB+nmutABf-nmutABi

                # Sortie précoce du protocole de boucle si i = 27
                if(i==27){
                    motherA = c(motherAnm,motherAhm)
                    motherB = c(motherBnm,motherBhm)
                    motherAB = c(motherABnm,motherABhm)
                }

                # Impression resultats à chaque generation
                    verboseA=sum(motherAnm)+sum(motherAhm)
                    verboseB=sum(motherBnm)+sum(motherBhm)
                    verboseAB=length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
                    print(paste("generation : ", i, "res : ",verboseA," ",verboseB,"",verboseAB," ", nmutA," ",nmutB," ",nmutAB))
            } 

############ generation de regime hypermutateur    
            else if (i>tgen && i<tgen+ngen) # traiter séparémenent les deux populations
            {
                # on garde ici les deux population séparée pour x générations
                # pour le calcul du nombre de mutants à cette generation
                nmutAi = 2*sum(motherAnm)+2*sum(motherAhm)
                nmutBi = 2*sum(motherBnm)+2*sum(motherBhm)
                nmutABi = 2*length(motherABnm[motherABnm>1])+2*length(motherABhm[motherABhm>1])

                # ici, on ne doit pas séparer ici les populations, mais continuer avec celles crées avant
                daughterAnm = sample(0:1,length(motherAnm),prob = c(1-x,x), replace = T)
                daughterBnm = sample(0:1,length(motherBnm),prob = c(1-y,y), replace = T)
                daughterAhm = sample(0:1,length(motherAhm),prob = c(1-x*m,x*m), replace = T)
                daughterBhm = sample(0:1,length(motherBhm),prob = c(1-y*m,y*m), replace = T)

                # trouver les doubles apparus dans les filles (chez les normaux et les hypermutateurs
                tempABnm=daughterAnm+daughterBnm
                tempABnm=replace(tempABnm, tempABnm>1,10)
                tempABnm=tempABnm+motherAnm+motherBnm # ajout des caractere maternel pour ne pas mettre des simultane dans les sequentiels
                candidaten=tempABnm[tempABnm==10]
                if(length(candidaten)>0)
                  { print(paste(length(candidaten)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }

                tempABhm=daughterAhm+daughterBhm
                tempABhm=replace(tempABhm, tempABhm>1,10)
                tempABhm=tempABhm+motherAhm+motherBhm # ajout pour ne pas mettre des simultane dans les sequentiels
                candidateh=tempABhm[tempABhm==10]
                if(length(candidateh)>0)
                  { print(paste(length(candidateh)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }

                # nettoyage des fichiers temp
                rm(tempABnm)
                rm(tempABhm)

                # transfert mother state to daughter state
                daughterAnm = motherAnm + daughterAnm
                daughterBnm = motherBnm + daughterBnm
                daughterAhm = motherAhm + daughterAhm
                daughterBhm = motherBhm + daughterBhm

                # junction of 2 "daughter" lists
                motherAnm = c(motherAnm, daughterAnm)
                motherBnm = c(motherBnm, daughterBnm) 
                motherAhm = c(motherAhm, daughterAhm)
                motherBhm = c(motherBhm, daughterBhm) 

                # correct for mutation in mother already mutated
                motherAnm=replace(motherAnm, motherAnm>0,1)
                motherBnm=replace(motherBnm, motherBnm>0,1)
                motherAhm=replace(motherAhm, motherAhm>0,1)
                motherBhm=replace(motherBhm, motherBhm>0,1)

                # table for dble mutant
                motherABnm=motherAnm+motherBnm
                motherABhm=motherAhm+motherBhm

                # nettoyage des fichiers temp
                rm(daughterAnm)        
                rm(daughterBnm)
                rm(daughterAhm)        
                rm(daughterBhm)

                # calcul du nombre de mutant à cette generation # ICI ON ADDITIONNE LES DEUX OU REMERGE LES DEUX POPULATIONS
                nmutAf = sum(motherAnm)+sum(motherAhm)
                nmutBf = sum(motherBnm)+sum(motherBhm)
                nmutABf = length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
                nmutA=nmutA+nmutAf-nmutAi
                nmutB=nmutB+nmutBf-nmutBi
                nmutAB=nmutAB+nmutABf-nmutABi

                # Sortie précoce du protocole de boucle si i=27 (si on a atteind le maximum de generation simulee)
                if(i==27){
                    motherA = c(motherAnm,motherAhm)
                    motherB = c(motherBnm,motherBhm)
                    motherAB = c(motherABnm,motherABhm)
                }

                # Impression resultats à chaque generation
                    verboseA=sum(motherAnm)+sum(motherAhm)
                    verboseB=sum(motherBnm)+sum(motherBhm)
                    verboseAB=length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
                    print(paste("generation : ", i, "res : ",verboseA," ",verboseB,"",verboseAB," ", nmutA," ",nmutB," ",nmutAB))
            }    

############ generation de sortie regime hypermutateur    
            else if (i==tgen+ngen) # sortir de la boucle et re-fusionner les populations
            {
                # pour le calcul du nombre de mutants à cette generation
                #nmutAi = 2*sum(motherA)
                #nmutBi = 2*sum(motherB)
                #nmutABi = 2*length(motherAB[motherAB>1])
                nmutAi = 2*sum(motherAnm)+2*sum(motherAhm)
                nmutBi = 2*sum(motherBnm)+2*sum(motherBhm)
                nmutABi = 2*length(motherABnm[motherABnm>1])+2*length(motherABhm[motherABhm>1])


                # ici, on ne doit pas séparer ici les populations, mais continuer avec celles crées avant
                daughterAnm = sample(0:1,length(motherAnm),prob = c(1-x,x), replace = T)
                daughterBnm = sample(0:1,length(motherBnm),prob = c(1-y,y), replace = T)
                daughterAhm = sample(0:1,length(motherAhm),prob = c(1-x*m,x*m), replace = T)
                daughterBhm = sample(0:1,length(motherBhm),prob = c(1-y*m,y*m), replace = T)

                # merge population (normal and transient hypermutator)
                daughterA=c(daughterAnm,daughterAhm) 
                daughterB=c(daughterBnm,daughterBhm)
                motherA=c(motherAnm,motherAhm) # mise à jour de la liste motherA qui prend en compte les 4 cycles d'hypermutateurs
                motherB=c(motherBnm,motherBhm) # mise à jour de la liste motherA qui prend en compte les 4 cycles d'hypermutateurs

                # search of double simu in new mutants
                tempAB=daughterA+daughterB
                tempAB=replace(tempAB, tempAB>1,10)
                tempAB=tempAB+motherA+motherB # ajout des caractere maternel pour ne pas mettre des simultane dans les sequentiels
                candidate=tempAB[tempAB==10]

                if(length(candidate)>0)
                  { print(paste(length(candidate)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }
         
                # nettoyage des fichiers temp
                rm(tempAB)

                # transfert mother state to daugther state
                daughterA = motherA + daughterA
                daughterB = motherB + daughterB
   
                # junction of 2 "daughter" lists
                motherA = c(motherA, daughterA)
                motherB = c(motherB, daughterB) 

                # correct for mutation in mother already mutated
                motherA=replace(motherA, motherA>0,1)
                motherB=replace(motherB, motherB>0,1)

                # table for dble mutant
                motherAB=motherA+motherB

                # nettoyage des fichiers temp
                rm(daughterA)        
                rm(daughterB)

                # calcul du nombre de mutant à cette generation
                nmutAf = sum(motherA)
                nmutBf = sum(motherB)
                nmutABf = length(motherAB[motherAB>1])
                nmutA=nmutA+nmutAf-nmutAi
                nmutB=nmutB+nmutBf-nmutBi
                nmutAB=nmutAB+nmutABf-nmutABi

                # Impression resultats à chaque generation
                print(paste("generation : ", i, "res : ",sum(motherA)," ",sum(motherB)," ",length(motherAB[motherAB>1])," ",nmutA," ",nmutB," ",nmutAB))
            }   

############ generations de regime normal    
            else
            {
                daughterA = sample(0:1,length(motherA),prob = c(1-x,x), replace = T)
                daughterB = sample(0:1,length(motherB),prob = c(1-y,y), replace = T)

                # pour le calcul du nombre de mutants à cette generation
                nmutAi = 2*sum(motherA)
                nmutBi = 2*sum(motherB)
                nmutABi = 2*length(motherAB[motherAB>1])

                # search of double simu in new mutants
                tempAB=daughterA+daughterB
                tempAB=replace(tempAB, tempAB>1,10)
                tempAB=tempAB+motherA+motherB # ajout des caractere maternel pour ne pas mettre des simultane dans les sequentiels
                candidate=tempAB[tempAB==10]

                if(length(candidate)>0)
                  { print(paste(length(candidate)," mutant simultanee a la generation :",as.character(i), " pour le cycle :",as.character(c)))
                  }
         
                # nettoyage des fichiers temp
                rm(tempAB)

                # transfert mother state to daugther state
                daughterA = motherA + daughterA
                daughterB = motherB + daughterB
    
                # junction of 2 "daughter" lists
                motherA = c(motherA, daughterA)
                motherB = c(motherB, daughterB) 

                # correct for mutation in mother already mutated
                motherA=replace(motherA, motherA>0,1)
                motherB=replace(motherB, motherB>0,1)

                # table for dble mutant
                motherAB=motherA+motherB

                # nettoyage des fichiers temp
                rm(daughterA)        
                rm(daughterB)

                # calcul du nombre de mutant à cette generation
                nmutAf = sum(motherA)
                nmutBf = sum(motherB)
                nmutABf = length(motherAB[motherAB>1])
                nmutA=nmutA+nmutAf-nmutAi
                nmutB=nmutB+nmutBf-nmutBi
                nmutAB=nmutAB+nmutABf-nmutABi

                # Impression resultats à chaque generation
                print(paste("generation : ", i, "res : ",sum(motherA)," ",sum(motherB)," ",length(motherAB[motherAB>1])," ",nmutA," ",nmutB," ",nmutAB))
            }
        }

    #------------------------------------
    # calcul nb colonies individuelles  |
    #------------------------------------

    resAB=motherAB[motherAB>1]
    A=sum(motherA)
    B=sum(motherB)
    C=length(resAB)
    reslist=c(A,B,C,nmutA,nmutB,nmutAB)

    # nettoyage des fichiers temp
    rm(motherA)        
    rm(motherB)
    rm(motherAB)

    #----------
    # output  |
    #----------
    return(reslist)
}


#################################
# Choix de l'espace de travail  #
#################################

# choix de dossier de travail
cat("\n")
myFolder=choose.dir(readline("Press <enter> to chose the current working folder\t"))
cat("\n")
setwd(myFolder)


#################################
# Boucle principale             #
#################################

continueWP = "y"
while (continueWP!= "n")
{

    ##################################
    # rentrer les donnees de depart  #
    ##################################

    # donner un nom à la simulation
    cat("\n")
    Output=as.character(readline("Choose a name for the result file :\t"))
    cat("\n")

    # definition taux 1
    r1=as.numeric(readline("Rate of mutation A :\t"))

    # definition taux 2
    r2=as.numeric(readline("Rate of mutation B :\t"))

    # definition nb de generation
    g=as.numeric(readline("Number of generation :\t"))

    #nb de repetitions
    n=as.numeric(readline("Number of repeat in the simulation :\t"))

    #taille population hypermut
    hyp=as.numeric(readline("Size of the transient population (nb cells) :\t"))

    #facteur augmentation taux
    mult=as.numeric(readline("Increase factor of hypermutator rates :\t"))

    #durée augmentation taux
    ngen=as.numeric(readline("Number of successive hypermutator generations :\t"))
    ngen=ngen-1 #retire la première génération qui n'est pas comprise dans le nombre de gen successive

    ##############################
    # 2 mutations independantes  #
    ##############################

    #-----------------------------------------
    # création des objets pour le résultats  |
    #-----------------------------------------

    res=matrix(nrow = n, ncol = 6)

    #--------------------------
    # Boucle de simulation    |
    #--------------------------

    sink(paste(Output,".independant.verbose.txt"), split = T)
    # summary taux
    print(paste("le taux choisi pour la mutation A est : ",r1))
    print(paste("le taux choisi pour la mutation B est : ",r2))
    print(paste("le nombre de generation est : ",g))
    print(paste(hyp," cellules hypermutatrice (x",mult,") sur ",ngen+1," generations")) 

    start=Sys.time()
    for (i in 1:n)
    {
        # verbose
        print(paste("simulation cycle : ",i," on ",n))

        # output results
        S=SIMUI(r1,r2,g,i,hyp,mult,tgen,ngen)
        res[i,1]=S[1]
        res[i,2]=S[2]
        res[i,3]=S[3] 
        res[i,4]=S[4]
        res[i,5]=S[5]
        res[i,6]=S[6]

        # calcul temp restant
        end=Sys.time()
        tdiff=difftime(end,start, units = "mins")                       
        print(paste("Il reste :", round(tdiff/i*(n-i),2)," mins"))
    }

    colnames(res)=c("A","B","AB","mutA","mutB","mutAB")
    write.table(res,paste(Output,".res_independant.txt"),row.names = FALSE, sep = "\t", quote = FALSE)
    sink()

    ####################################
    # Relancer une autre simulation    #
    ####################################

    cat("\n")
    Hello=as.character(readline("Do you want to perform another simulation (y or n):\t"))
    cat("\n")

    #petite sécurité en cas de frappe non valide
    while ((Hello != "y") && (Hello != "n"))
        {
            cat("Please choose yes <y> or no <n>\n")
            cat("\n")
            Hello=as.character(readline("Do you want to perform another simulation (y or n):\t"))
            cat("\n") 
        }

    continueWP=Hello
}

cat("Have a nice day :-)\n")
cat("\n")
