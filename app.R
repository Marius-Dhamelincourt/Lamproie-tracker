
#Package loading
library(shiny)
library(shinycssloaders)
library(MASS)
library(EasyABC)
library(abc)
library(abctools)
#install.packages("countreg", repos="http://R-Forge.R-project.org")
#library(countreg)# for zero truncated negative binomial distributions
library(msm)# for zero truncated normal distribution
library(fGarch)#Compute skewness normal distribution
library(readxl)#to read excel files
library(fitdistrplus)#Distributions fitting
library(DT)
library(reshape2)#to use melt function
library(ggthemes)#to use theme stata in ggplot graphs
library(shinydashboard)# to make boxes
library(shinybusy)# To add loading icon

# Function to manage different dataframe length in plot (due to sex-ratio different from 0.5)
cbind_dif <- function(x = list()){
  # Find max length
  max_length <- max(unlist(lapply(x, length)))
  
  # Set length of each vector as
  res <- lapply(x, function(x){
    length(x) <- max_length
    return(x)
  })
  
  return(as.data.frame(res))
}

#Function to calculate zero-truncated binomial distribution. Initially provided within the "countreg"
#package in R-Forge repository, but provided here because of the non-availability of R-Forge packages for Shiny apps
dztnbinom <- function(x, mu, theta, size, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) - pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[mu <= 0] <- 0
  if(log) rval else exp(rval)
}

qztnbinom <- function(p, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dnbinom(0, mu = mu, size = theta)
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztnbinom(1, mu = mu, theta = theta, log = log.p)] <- 1
  rval
}

rztnbinom <- function(n, mu, theta, size) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qztnbinom(runif(n), mu = mu, theta = theta)
}

#UI
##################################################################################################################################################################
##################################################################################################################################################################
#UI

ui<-fluidPage(navbarPage(inverse=TRUE,"Lamproie tracker",
           tabPanel("Comment utiliser cette application ?",
           fluidRow(
             p(style="text-align: justify;","Cette application a été développée afin de pouvoir déterminer une population de lamproies marines
                présente sur un site au cours d'une saison de reproduction. Pour ce faire, l'utilisateur doit fournir
                un jeu de données au format .csv dans l'onglet 'Chargement des données' avec en ligne les différents 
                jours de suivi (une ligne = un jour) et en colonnes le nombre de nids actifs recensés chaque jour."), 
             p(style="text-align: justify;","Un nid actif est un nid sur lequel au moins un individu, mâle ou femelle, est présent. Lors d'un comptage,
                l'utilisateur recense les nids indépendamment de ce qu'il a observé les jours précédents.
                Ainsi, si un nid était déjà observé la veille, il est de nouveau compté s'il est occupé. 
                L'application simule une saison de reproduction en tenant compte de différents paramètres 
                fixés au préalable (Voir onglet 'Paramètres'). Ceux-ci ont été obtenus au cours d'une 
                saison de reproduction lors de laquelle des individus ont été marqués et suivis par un 
                protocole de Capture-Marquage-Recapture couplé à des observations comportementales. 
                L'utilisateur peut visualiser ces paramètres et éventuellement les modifier s'il le souhaite.
                Cependant, il est conseillé de ne pas les changer sans raison, ceux-ci étant directement utilisés par le modèle."),
             p(style="text-align: justify;"," A la suite du chargement des données, l'analyse peut directement être lancée dans l'onglet 'Lancement de l'analyse'.
                Celle-ci peut durer plusieurs minutes et dépend des paramètres de vos données. Lorsque celle-ci est terminée,
                un graphique de distribution de la valeur de population estimée est affiché, ainsi qu'un tableau indiquant 
                la médiane et les quantiles à 2,5 et 97,5% (soit 95% de confiance) de la distribution estimée.")),
           column(6,titlePanel(div(img(src = "paysage.png", width = "100%", class = "bg"),))),
           column(6,titlePanel(div(img(src = "lamproie_substrat.png", width = "100%", class = "bg"),)))),
           tabPanel("Chargement des données",
                    fluidRow(
                      column(12,wellPanel(h3("Import"),
                                         fileInput("file1", "Choisissez le fichier csv",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                                         checkboxInput("header", "Mes données ont des en-têtes", TRUE),
                                         # Input: Select separator ----
                                         radioButtons("sep", "Séparateur",
                                                      choices = c(Virgule = ",","Point virgule" = ";",Espace = "\t"),
                                                      selected = ";"),
                                         # Input: Select decimal separator ----
                                         radioButtons("dec", "Marqueur de décimale",
                                                      choices = c(Virgule = ",","Point" = "."),
                                                      selected = "."),
                                         # Input: Select quotes ----
                                         radioButtons("quote", "Guillemets",
                                                      choices = c("Pas de guillemets" = "","Guillemets doubles" = '"',"Guillemets simples" = "'"),
                                                      selected = '"')
                                         
                                         
                                         
                      )
                      )
                      
                    ),
                fluidRow(
                  
                column(12, wellPanel(h3("Valeurs a priori"),
                                     numericInput("min.prior","Nombre minimal de lamproies, a priori",value=10),
                                     numericInput("max.prior","Nombre maximal de lamproies, a priori",value=250)))),
                fluidRow(tableOutput("contents")) 
           ),
           tabPanel("Paramètres",
                    fluidRow(
                      column(12,wellPanel(h4("Les paramètres choisis par défaut ont été obtenus au cours d'une étude réalisée sur la Nive en 2019 par suivi de 115 individus sur un site au cours d'une saison de reproduction (Dhamelincourt et al., 2021). Pour en savoir plus consulter les liens suivants : "),
                    
                    fluidPage(uiOutput("tab")),
                    fluidPage(uiOutput("tab1"))))),
                    
                    fluidRow(wellPanel(
                        #numericInput("n_individuals","Nombre d'individus (uniquement utilisé pour visualiser les paramètres)",value=200,step=1),
                        numericInput("sex.ratio","Sex-ratio",value=0.5,step=0.05))),
                   
                    fluidRow(
                      column(12, h3("Nombre d'individus par nid"), h4("Ce paramètre indique la distribution attendue du nombre de mâles et du nombre de femelles occupant un nid donné. Une distribution binomiale négative tronquée en 0 a été appliquée."))),
                       
                    
                    fluidRow(column(12, 
                      box(wellPanel( 
                        numericInput("mu.number.ind.nest.male","Nombre de mâles par nid : moyenne",value=1.25,step=0.01),
                        numericInput("size.number.ind.nest.male","Nombre de mâles par nid : paramètre de dispersion",value=407.28,step=0.01),
                        numericInput("mu.number.ind.nest.female","Nombre de femelles par nid : moyenne",value=1.15,step=0.01),
                        numericInput("size.number.ind.nest.female","Nombre de femelles par nid : paramètre de dispersion",value=4.55,step=0.01))),
                     box(plotOutput("number_individuals")))),
                    
                    fluidRow(
                      column(12, h3("Nombre de nids par individu"),h4("Ce paramètre indique le nombre de nids que fréquente un individu au cours de sa saison de reproduction. Une distribution binomiale négative tronquée en 0 a été appliquée."))),
                    
                    
                    fluidRow(column(12, 
                      box(wellPanel( 
                          numericInput("mu.number.nest.ind.male","Nombre de nids par mâle : moyenne",value=2.38,step=0.01), 
                          numericInput("size.number.nest.ind.male","Nombre de nids par mâle : paramètre de dispersion",value=21.13,step=0.01),
                          numericInput("mu.number.nest.ind.female","Nombre de nids par femelle : moyenne",value=1.70,step=0.01),
                          numericInput("size.number.nest.ind.female","Nombre de nids par femelle : paramètre de dispersion",value=305.15,step=0.01))),
                      box(plotOutput("number_nests")))),
                    
                    fluidRow(
                      column(12, h3("Délai d'arrivée sur le site"),h4("Le délai d'arrivée correspond à l'écart en jours entre le début de l'activité de reproduction sur le site et l'arrivée de l'individu. Une loi normale asymétrique tronquée en 0 a été utilisée."))),
                      
                    fluidRow(column(12, 
                                    box(wellPanel( 
                        numericInput("mu.delay.male","Moyenne du délai d'arrivée sur le site pour les mâles",value=27.25,step=0.01),
                        numericInput("sd.delay.male","Ecart-type du délai d'arrivée sur le site pour les mâles",value=11.04,step=0.01),
                        numericInput("mu.delay.female","Moyenne du délai d'arrivée sur le site pour les femelles",value=30.45,step=0.01),
                        numericInput("sd.delay.female","Ecart-type du délai d'arrivée sur le site pour les femelles",value=10.31,step=0.01))),
                    box(plotOutput("delay")))),
                    
                    fluidRow(
                      column(12, h3("Temps de résidence"),h4("Le temps de résidence est la durée que chaque individu passe sur le site de reproduction après son délai d'arrivée, et au cours duquel il peut se reproduire. Une loi normale tronquée en 0 est ici utilisée."))),
                    
                    fluidRow(column(12, 
                                    box(wellPanel( 
                        numericInput("mean.residence.time.male","Moyenne du temps de résidence des mâles",value=8.33,step=0.01),
                        numericInput("sd.residence.time.male","Ecart-type du temps de résidence des mâles",value=1.02,step=0.01),
                        numericInput("mean.residence.time.female","Moyenne du temps de résidence des femelles",value=3.57,step=0.01),
                        numericInput("sd.residence.time.female","Ecart-type du temps de résidence des femelles",value=1.04,step=0.01))),
                        box(plotOutput("residence")))),
                    
                    fluidRow(
                      column(12, h3("Temps d'occupation d'un nid"),h4("Ce temps d'occupation correspond à la durée de présence d'un individu sur un nid. Une loi normale tronquée en 0 est utilisée."))),
                    
                    fluidRow(column(12, 
                                    box(wellPanel(
                        numericInput("mean.duration.nest","Temps d'occupation moyen ",value=1,step=0.01),
                        numericInput("sd.duration.nest","Ecart-type du temps d'occupation moyen ",value=0.5,step=0.01))),
                        box(plotOutput("duration")))),
                      ),
           #tabPanel("Lancement de l'analyse"),
           tabPanel("Lancement de l'analyse",
                    column(12,wellPanel(h4("Merci de vérifier que vous avez bien chargé les données au préalable dans l'onglet 'Chargement des données'. Une fois lancée, l'analyse peut prendre plusieurs minutes. Merci de patienter."))),
                    actionButton("launch", "Lancer l'analyse"),
                    fluidRow(
                      plotOutput("post.plot"),
                      tableOutput("abc.output"),
                      fluidRow(column(6,uiOutput("duree1"))),
                      )),
           tabPanel("Contacts et remerciements",
                    fluidRow(
                      column(12,wellPanel(h3("Contact des développeurs :"),
                             h4("Marius Dhamelincourt (Doctorant INRAE/UPPA/UPV) : marius.dhamelincourt@gmail.com"),
                             h4("Cédric Tentelier (Maître de Conférences INRAE/UPPA) : cedric.tentelier@univ-pau.fr")))),
                    fluidRow(h3("Partenaires et financeurs"),
                             h4(em("(Cliquez sur le logo pour en savoir plus)")),
                      column(4,tabItem("icratio",
                                      fluidRow(
                                        tags$a(img(src="ecobiop.png", height = 140, width = 250), href="https://ecobiop.com/", target = "_blank"))),
                               tabItem("icratio",
                                            fluidRow(
                                              tags$a(img(src="inrae.jpg", height = 140, width = 290), href="https://www.inrae.fr/", target = "_blank")))),
                      
                      column(4,tabItem("icratio",
                                      fluidRow(
                                         tags$a(img(src="uppa.jpg", height = 140, width = 210), href="https://www.univ-pau.fr/fr/index.html", target = "_blank"))),
                                        
                               tabItem("icratio",
                                                fluidRow(
                                                  tags$a(img(src="upv.png", height = 140, width = 250), href="https://www.ehu.eus/es/home", target = "_blank")))),
                      column(4,tabItem("icratio",
                                       fluidRow(
                                           tags$a(img(src="ofb.jpg", height = 140, width = 210), href="https://www.ofb.gouv.fr/", target = "_blank"))),
                                 
                                    tabItem("icratio",
                                         fluidRow(
                                           tags$a(img(src="agrocampus.png", height = 140, width = 250), href="http://formationcontinue.agrocampus-ouest.fr/infoglueDeliverLive/", target = "_blank"))))
                             
                       ) 
                    
           
        )
))


#SERVER
##################################################################################################################################################################
##################################################################################################################################################################
#SERVER

server<-function(input, output, session) {

##Hyperlink to the article

url <- a("Lien vers l'article (Journal of Fish Biology)", href="https://onlinelibrary.wiley.com/doi/10.1111/jfb.14601", target = "_blank")
output$tab <- renderUI({
    tagList("", url)
})

url1 <- a("Lien vers le site de dépôt (HAL INRAe)", href="https://hal.inrae.fr/hal-03036001", target = "_blank")
output$tab1 <- renderUI({
  tagList("", url1)
})


  
##Data loading
  
data_lamprey<-reactive(read.csv(input$file1$datapath,
                          header = input$header,
                          sep = input$sep,
                          dec = input$dec,
                          quote = input$quote))
  
output$contents <- renderTable({
    
# input$file1 will be NULL initially. After the user selects
# and uploads a file, head of that data file by default,
# or all rows if selected, will be shown.
req(input$file1)
    
# when reading semicolon separated files,
# having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- read.csv(input$file1$datapath,
                       header = input$header,
                       sep = input$sep,
                       dec = input$dec,
                       quote = input$quote)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    return(df)
    
  })

###Plots of parameter distribution

#Number of individuals in each nest

data_to_add1 <- reactive({
  nb_individuals_m<- rztnbinom(1000,size=input$size.number.ind.nest.male,mu=input$mu.number.ind.nest.male)
  nb_individuals_f<- rztnbinom(1000,size=input$size.number.ind.nest.female,mu=input$mu.number.ind.nest.female)
  data_individuals<-melt(cbind_dif(list(nb_individuals_m=nb_individuals_m, nb_individuals_f=nb_individuals_f)))
})

output$number_individuals <- renderPlot({
  data_individuals<-data_to_add1()
  ggplot(data_individuals, aes(x = value, group = variable)) +
    geom_bar(aes(fill = variable), alpha = 0.8, position = 'identity')+
    theme_stata(scheme = "s1mono")+
    theme(axis.title.x=element_text(size=12,face="bold",vjust=-2),
          axis.title.y=element_text(size=12,face="bold",vjust=5),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"),
          strip.text.x = element_blank())+
    xlab(label="Nombre d'individus par nid")+
    ylab(label="Comptage")+
    scale_fill_manual(values=c("#0072B2","#CC79A7"), labels= c("Mâles", "Femelles"))+ 
    labs(fill = "Sexe")+
    scale_x_continuous(breaks = scales::breaks_width(1))+
    facet_wrap(~variable,  ncol=1)
})

#Number of nests for each individual

data_to_add<- reactive({
  nb_nests_m <-rztnbinom(1000,size=input$size.number.nest.ind.male,mu=input$mu.number.nest.ind.male)
  nb_nests_f <-rztnbinom(1000,size=input$size.number.nest.ind.female,mu=input$mu.number.nest.ind.female)
  data_nests <- melt(cbind_dif(list(nb_nests_m=nb_nests_m, nb_nests_f=nb_nests_f)))
})

output$number_nests <- renderPlot({
data_nests<-data_to_add()
ggplot(data_nests, aes(x = value, group = variable)) +
  geom_bar(aes(fill = variable), alpha = 0.8, position = 'identity')+
  theme_stata(scheme = "s1mono")+
  theme(axis.title.x=element_text(size=12,face="bold",vjust=-2),
        axis.title.y=element_text(size=12,face="bold",vjust=5),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        strip.text.x = element_blank())+
  xlab(label="Nombre de nids par individu")+
  ylab(label="Comptage")+
  scale_fill_manual(values=c("#0072B2","#CC79A7"), labels= c("Mâles", "Femelles"))+ 
  labs(fill = "Sexe")+
  scale_x_continuous(breaks = scales::breaks_width(1))+
  facet_wrap(~variable,  ncol=1)
})





#Delay of arrival on spawning ground

data_to_add2 <- reactive({
delay.m<-round(rsnorm(2000*10, mean = input$mu.delay.male, sd = input$sd.delay.male, xi = -0.94))#Delay of arrival on the spawning ground for males
delay.m<-delay.m[which(delay.m > 0)]#we select positive values
delay.m<- sample(delay.m, 1000, replace=TRUE)# we draw Nmales positive values within the normal skewed distribution

delay.f<-round(rsnorm(2000*10, mean = input$mu.delay.female, sd = input$sd.delay.female, xi = -1.40))#Delay of arrival on the spawning ground for females
delay.f<-delay.f[which(delay.f > 0)]#we select positive values
delay.f<- sample(delay.f, 1000, replace=TRUE)# we draw Nmales positive values within the normal skewed distribution
data_delay<-melt(cbind_dif(list(delay.m=delay.m, delay.f=delay.f)))
})

output$delay <- renderPlot({
  data_delay<-data_to_add2()
  ggplot(data_delay, aes(x = value, group = variable)) +
    geom_bar(aes(fill = variable), alpha = 0.8, position = 'identity')+
    theme_stata(scheme = "s1mono")+
    theme(axis.title.x=element_text(size=12,face="bold",vjust=-2),
          axis.title.y=element_text(size=12,face="bold",vjust=5),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"),
          strip.text.x = element_blank())+
    xlab(label="Nombre de jours de retard")+
    ylab(label="Comptage")+
    scale_fill_manual(values=c("#0072B2","#CC79A7"), labels= c("Mâles", "Femelles"))+ 
    labs(fill = "Sexe")+
    scale_x_continuous(breaks = scales::breaks_width(10))+
    facet_wrap(~variable,  ncol=1)
})

#Residence time

 data_to_add3 <- reactive({
 residence.m<-round(rnorm(1000,mean=input$mean.residence.time.male,sd=input$sd.residence.time.male)) 
 residence.f<-round(rnorm(1000,mean=input$mean.residence.time.female,sd=input$sd.residence.time.female))
 data_residence<-melt(cbind_dif(list(residence.m=residence.m, residence.f=residence.f)))
 data_residence<-data_residence[data_residence$value >= 1,]
 })
 
 output$residence <- renderPlot({
   data_residence<-data_to_add3()
   ggplot(data_residence, aes(x = value, group = variable)) +
     geom_bar(aes(fill = variable), alpha = 0.8, position = 'identity')+
     theme_stata(scheme = "s1mono")+
     theme(axis.title.x=element_text(size=12,face="bold",vjust=-2),
           axis.title.y=element_text(size=12,face="bold",vjust=5),
           axis.text.x =element_text(size=10,face="bold"),
           axis.text.y =element_text(size=10,face="bold"),
           strip.text.x = element_blank())+
     xlab(label="Temps de résidence")+
     ylab(label="Comptage")+
     scale_fill_manual(values=c("#0072B2","#CC79A7"), labels= c("Mâles", "Femelles"))+ 
     labs(fill = "Sexe")+
     scale_x_continuous(breaks = scales::breaks_width(1))+
     facet_wrap(~variable,  ncol=1)
 })

#Duration on nest

data_to_add4 <- reactive({
duration<-round(rtnorm(2000,mean=input$mean.duration.nest,sd=input$sd.duration.nest, lower=1))
data_duration<-data.frame(duration)
})

output$duration <- renderPlot({
  data_duration<-data_to_add4()
  ggplot(data_duration, aes(x = duration)) +
    geom_bar(fill = "#009E73", alpha = 0.8, position = 'identity')+
    theme_stata(scheme = "s1mono")+
    theme(axis.title.x=element_text(size=12,face="bold",vjust=-2),
          axis.title.y=element_text(size=12,face="bold",vjust=5),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"))+
    xlab(label="Délai entre chaque construction de nid")+
    ylab(label="Comptage")+
    scale_x_continuous(breaks = scales::breaks_width(1))
})

##Sea lamprey model
lamprey_spawning<-function(N,
                           sex.ratio= input$sex.ratio,
                           size.number.ind.nest.male=input$size.number.ind.nest.male, mu.number.ind.nest.male=input$mu.number.ind.nest.male,  #Number of males per nest
                           size.number.ind.nest.female=input$size.number.ind.nest.female, mu.number.ind.nest.female=input$mu.number.ind.nest.female, #Number of females per nest
                           
                           size.number.nest.ind.male=input$size.number.nest.ind.male , mu.number.nest.ind.male=input$mu.number.nest.ind.male, #Number of nests per male
                           size.number.nest.ind.female=input$size.number.nest.ind.female ,mu.number.nest.ind.female=input$mu.number.nest.ind.female, #Number of nests per female
                           
                           mu.delay.male=input$mu.delay.male, sd.delay.male=input$sd.delay.male,  #Delay of arrival on the spawning ground for males: mu.delay.male=27.25, sd.delay.male=11.04
                           mu.delay.female=input$mu.delay.female, sd.delay.female=input$sd.delay.female,  #Delay of arrival on the spawning ground for females: mu.delay.female=30.45, sd.delay.female=10.31
                           
                           mean.residence.time.male=input$mean.residence.time.male, sd.residence.time.male=input$sd.residence.time.male,  #Residence time for males
                           mean.residence.time.female=input$mean.residence.time.female, sd.residence.time.female=input$sd.residence.time.female,  #Residence time for females
                           
                           mean.duration.nest= input$mean.duration.nest, sd.duration.nest= input$sd.duration.nest) { #Time between nests calculated only with nests > 1 individual.
  
  
  time<-1:100 #we can consider a 100 days spawning season (to be large)
  
  
  
  ######## Number of individuals for each sex ########
  Nmales=round(N*sex.ratio)
  Nfemales=round(N*(1-sex.ratio))
  
  ######## Delay males and females ########
  
  delay.male<-round(rsnorm(N*10, mean = mu.delay.male, sd = sd.delay.male, xi = -0.94))#Delay of arrival on the spawning ground for males
  delay.male<-delay.male[which(delay.male > 0)]#we select positive values
  delay.male<- sample(delay.male, Nmales, replace=TRUE)# we draw Nmales positive values within the normal skewed distribution
  
  delay.female<-round(rsnorm(N*10, mean = mu.delay.female, sd = sd.delay.female, xi = -1.40))#Delay of arrival on the spawning ground for females
  delay.female<-delay.female[which(delay.female > 0)]#we select positive values
  delay.female<- sample(delay.female, Nfemales, replace=TRUE)# we draw Nmales positive values within the normal skewed distribution
  
  
  ######## Residence time ########  
  residence.male<-round(rtnorm(Nmales,mean=mean.residence.time.male,sd=sd.residence.time.male, lower=1)) 
  residence.female<-round(rtnorm(Nfemales,mean=mean.residence.time.female,sd=sd.residence.time.female, lower=1)) # Residence time
  
  ######## Duration on the nest ########
  duration.nest<-round(rtnorm(Nmales,mean=mean.duration.nest,sd=sd.duration.nest, lower=1)) 
  
  ########Simulate a number of nests for each individual########
  Nb.nests.female<-rztnbinom(Nfemales,size=size.number.nest.ind.female,mu=mu.number.nest.ind.female)
  Nb.nests.male<-rztnbinom(Nmales,size=size.number.nest.ind.male,mu=mu.number.nest.ind.male)
  
  ########Simulate a number of individuals per nest#######
  Nb.female.nest<- rztnbinom(Nfemales,size=size.number.ind.nest.female,mu=mu.number.ind.nest.female)
  Nb.male.nest<- rztnbinom(Nmales,size=size.number.nest.ind.male,mu=mu.number.nest.ind.male)
  
  ########Creation of empty arrays for males and females nest count and total nest count########
  
  Nb.of.nests.male<-array(dim=c(length(1:Nmales)),0)
  Nb.of.nests.female<-array(dim=c(length(1:Nfemales)),0)
  Nb.of.nests.total<-array(dim=c(length(time)),0)
  Nb.of.nests.one.male<-array(dim=c(length(time)),0)# nests dug by only one male
  Nb.of.nests.one.female<-array(dim=c(length(time)),0)# nests dug by only one female
  
  ########Matrix of males and females active during day i#######
  
  spawning.day<-array(dim=c(Nmales,Nfemales,length(time)),data=0)
  
  
  
  for(i in 1:length(time)){
    for(j in sample(c(1 : Nmales))){ #We randomize the iteration order to ensure the individuals with high index are not always considered after the others
      for(k in sample(c(1 : Nfemales))){
        
        if(delay.female[k]<i && delay.male[j]<i && (delay.female[k] + residence.female[k])>i && (delay.male[j] + residence.male[j]) > i) 
        {spawning.day[j,k,i]<-1}
        
      }
    }
  }
  
  spawning.day.initial<- spawning.day #Used to check if the behaviour of the loop is fine
  
  ########Arrange individuals in nest respecting number of nests per individual#######
  
  for(i in 1:length(time)){ 
    
    random_Nb.female.nest<-sample(Nb.female.nest,length(Nb.female.nest),replace=FALSE)#We choose randomly one number of partners of opposite sex each day
    random_Nb.male.nest<-sample(Nb.male.nest,length(Nb.male.nest),replace=FALSE)#We choose randomly one number of partners of opposite sex each day
    
    for(j in 1 : Nmales){ #We randomize the iteration order to ensure the individuals with high index are not always considered after the others
      for(k in 1 : Nfemales){
        
        if (spawning.day[j,k,i]==1 && 
            matrixStats::count(spawning.day[j,,i], value = 2) < random_Nb.female.nest[k]&& # each day an individual cannot make nest with more individuals than the limit defined
            matrixStats::count(spawning.day[,k,i], value = 2) < random_Nb.male.nest[j]&&
            (i+duration.nest[j]) <= length(time))# to avoid a script error and nests built after the end of the season
          
        {spawning.day[j,k,i:(i+duration.nest[j])]<- 2 # a duration time is added to take into account the time of occupation for a nest.
        
        } 
        
      }
    }
  }   
  
  
  
  
  ######### Return the total number of nests and the number of nests for each individual ########
  
  for(i in 1:length(time)){
    for(j in 1 : Nmales){ 
      
      
      for(k in 1 : Nfemales){
        
        if (spawning.day[j,k,i]== 2 && "2" %in% spawning.day[(1:j-1),k,i] == FALSE && "2" %in% spawning.day[j,(1:k-1),i] == FALSE) #&& "2" %in% spawning.day[j,k,(i-1)] == FALSE )
          
        {Nb.of.nests.total[i] <- Nb.of.nests.total[i]+1}
        
         if (spawning.day[j,k,i]== 2 && "2" %in% spawning.day[j,k,i-1] == FALSE && "2" %in% spawning.day[j,(1:k-1),i] == FALSE)
         {Nb.of.nests.male[j] <- Nb.of.nests.male[j]+1}
         
         if (spawning.day[j,k,i]== 2 && "2" %in% spawning.day[j,k,i-1] == FALSE && "2" %in% spawning.day[(1:j-1),k,i] == FALSE)
         {Nb.of.nests.female[k] <- Nb.of.nests.female[k]+1}
        
        
        
      }
    }
  }
  
  for(i in 1:length(time)){
    for(j in 1 : Nmales){ 
      
      if ("1" %in% spawning.day[j,,i] == TRUE && 
          Nb.of.nests.male[j] < Nb.nests.male[j]) 
        #If a male have no partner during one day he digs a nest alone if he don't exceed a realistic number of nests
      {Nb.of.nests.one.male[i] <- Nb.of.nests.one.male[i]+1
      Nb.of.nests.male[j] <- Nb.of.nests.male[j] +1}
      
    }
  }
  
  for(i in 1:length(time)){
    for(k in 1 : Nfemales){ 
      
      if ("1" %in% spawning.day[,k,i] == TRUE && 
          Nb.of.nests.female[k] < Nb.nests.female[k]) 
        #If a female have no partner during one day he digs a nest alone if he don't exceed a realistic number of nests
      {Nb.of.nests.one.female[i] <- Nb.of.nests.one.female[i]+1
      Nb.of.nests.female[k] <- Nb.of.nests.female[k] +1}
      
    }
  }
  
  
  total_count <- Nb.of.nests.total + Nb.of.nests.one.male + Nb.of.nests.one.female
  
  #Now we must return the statistics for the total nest count from the first nest occurrence to the last (otherwise statistics are false)
  #If no nest is created it returns the full vector (to avoid problems with low populations)
  
  
  if(sum(total_count)!=0){
    indx <- which(total_count !=0)
    total_count<-total_count[indx[1L]:indx[length(indx)]]#we need to remove 0 before and after the spawning season, to calculate non-biased summary statistics
    
    
  } else{
    total_count# just in case we have only 0, the total dataset with 0 is returned
  }
  
  
  return(c( max(total_count), median(total_count), mean(total_count), quantile(total_count, 0.25),
            quantile(total_count, 0.75))) # The function returns a summary statictics vector
  
}


######## Function to use for ABC algorithm ########

lamprey_spawning.4ABC<-function(n.individuals){
  lamprey_spawning(N=as.integer(n.individuals))
}

#ABC model

observeEvent(input$launch,{
  show_modal_spinner(
    spin = "flower",
    color = "forestgreen",
    text = "Estimation en cours, merci de patienter..."
  )
  
# Add a prior for population size 
prior<-list(c("unif",input$min.prior,input$max.prior))#Distribution corresponding to the number of individuals REALLY observed (corresponding to the number of occupied nests observed, not the real number of nests)
# Indicate summary statistics 
summary.stats<-c(max(data_lamprey()[,2], na.rm=TRUE),median(data_lamprey()[,2], na.rm=TRUE),mean(data_lamprey()[,2], na.rm=TRUE), quantile(data_lamprey()[,2], 0.25, na.rm=TRUE), quantile(data_lamprey()[,2], 0.75, na.rm=TRUE) )

duree<-system.time({  
abc.seq<-ABC_sequential(method="Lenormand",model=lamprey_spawning.4ABC,prior=prior,nb_simul=50,summary_stat_target = summary.stats, progress_bar=TRUE)
})

#We calculate the duration of the calculation
duree1<-as.numeric(duree)[3]
output$duree1<-renderText({ 
  paste("Durée totale de calcul :"," ",duree1," ","s"," ","/"," ",duree1/60," ","min",sep="")
})

remove_modal_spinner()

d<-density(abc.seq$param,weights=abc.seq$weights)
x<-data.frame(t(wtd.quantile(x=abc.seq$param, weights=abc.seq$weights*1000, probs=c(0.025,0.5,0.975))))
colnames(x)<-c('Quantile 0.025','Mediane','Quantile 0.975')

output$post.plot<-renderPlot({
  plot(c(input$min.prior,input$max.prior),dunif(x=c(input$min.prior,input$max.prior),min=input$min.prior,max=input$max.prior)
       ,type="l",lwd=2,lty=3,ylim=c(0,max(d[['y']])),main="Distribution",xlab="Nombre d'individus",ylab="Probabilité")
  lines(d,type="l",col=2,lwd=2)
  abline(v=x,col="blue",lwd=c(1,2,1),lty=2)
  legend(x=ifelse(x$Mediane<mean(c(input$min.prior,input$max.prior)),"topright","topleft"),legend=c("a priori","a posteriori","Mediane","Quantiles 0.025 & 0.975"),col=c("black","red",rep("blue",3)),lwd=c(rep(2,3),1),lty=c(3,1,2,2))
})
output$abc.output<-renderTable({
  x
})

  

  


})



}

#APP LAUNCHING
##################################################################################################################################################################
##################################################################################################################################################################
#APP LAUNCHING

shinyApp(ui=ui, server=server)



