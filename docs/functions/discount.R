cat("discount - Função p/ o calculo de correçãode umidade e impureza  grãos
      Wilhan Valasco , contato: (64) - 992862039")
discount <- function(treat,prod,Ii,If,Ui,Uf){
      Trat <- factor(treat)
      prod <- as.numeric(prod)
      Prod <- prod
      # Prod (kg)
      # Impurza final e inicial (%)
      # Umidade Inicial e Final(%)
       QI <- ((Ii - If)*100/(100-If))
       QU <- ((Ui - Uf)*100/(100-Uf))
       DE <- 100 - ((100-QI)*(100-QU)/100)
       Discount <- ((Prod*DE)/100)
       Result <-  (Prod - Discount)
       Treatment <- treat
       Treatment <- factor(Treatment)
       Prod.uncorrected <- prod
       Prod.corrected <- Result
       return(data.frame(Treatment ,Prod.uncorrected,Prod.corrected))
  }


















