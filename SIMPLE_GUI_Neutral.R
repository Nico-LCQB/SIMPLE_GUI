cat("\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("%                                                                                       %\n")
cat("%       SIMPLE : SImulation de Mutations doubles dans des Populations de LEvures        %\n")
cat("%                                                                                       %\n")
cat("%                     nicolas Agier - V7.10 (Lang model-09-12-2022)                     %\n")
cat("%                                                                                       %\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("\n")
cat("\n")

# dans ce modele, seule l'une des deux filles issues de chaque division peut porter une mutation.
# Lang model
# dans cette version, la population hypermutatrice transitoire est simulée
# Comme inclus dans la simulation de deux mutations independantes, la simulation de 2 simultanées est inactivee
# Cette version donne le nombre de mutations ayant eu lieu dans chaque realisation

###########################
# Fonction de simulation  #
###########################

######################################
# fonction 2 mutations independantes #
######################################

SIMUI <- function(x,y,z,c,h,m,tgen) # x=r1, y=r2, z=g, c = cycle, h = size transient pop, m = increase of rates, tgen = generation of transiant hyperm
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
        # print(tgen)
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
            #print(i)
            if(i==tgen)
            {

                # pour le calcul du nombre de mutants à cette generation
                nmutAi = 2*sum(motherA)
                nmutBi = 2*sum(motherB)
                nmutABi = 2*length(motherAB[motherAB>1])

                # separer les populations hyper mutatrices et les autres
                daughterAnm = sample(0:1,length(motherA)-h,prob = c(1-x,x), replace = T)
                daughterBnm = sample(0:1,length(motherB)-h,prob = c(1-y,y), replace = T)
                daughterAhm = sample(0:1,h,prob = c(1-x*m,x*m), replace = T)
                daughterBhm = sample(0:1,h,prob = c(1-y*m,y*m), replace = T)

                # merge population (normal and transient hypermutator)
                daughterA=c(daughterAnm,daughterAhm) 
                daughterB=c(daughterBnm,daughterBhm)

                # search of double simu in new mutants
                tempAB=daughterA+daughterB
                tempAB=replace(tempAB, tempAB>1,10)
                tempAB=tempAB+motherA+motherB
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
            }
        
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
                tempAB=tempAB+motherA+motherB
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

    start=Sys.time()
    for (i in 1:n)
    {
        # verbose
        print(paste("simulation cycle : ",i," on ",n))

        # output results
        S=SIMUI(r1,r2,g,i,hyp,mult,tgen)
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
