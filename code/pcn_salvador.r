library(rsalvador)

setwd("~/SYNC_RPM/SYNC_Projects/Dominatrix/code")

#dataPath<-'/Users/ESB/ESB_DATA/DOMINATRIX_data/runs'
dataPath<-'../data/runs'

dirName<-'neutralModel_mutRate1e-08_Levels25'
nH=20
Hs<- seq(1/nH,1,length=nH)
PCNs<-c(1, 2, 20,200)


outFile<-paste(c(dataPath,dirName,"MLmut.csv"),collapse="/")
if (file.exists(outFile)) 
  file.remove(outFile)

for (PCN in PCNs){
	for (H in Hs){

		fileName<-paste(c('sim_PCN',PCN,'_H',H*100,'e-2.txt'), collapse="")
    print(fileName)

		dataFile <- c(dataPath, dirName,'data', fileName)
		if (file.exists(dataFile)) {
			y=import.text.data(paste(dataFile, collapse="/"))
			m=newton.LD(y)
			ci=confint.LD(y)

			data_row<-cbind(PCN=PCN,H=H,m=m,ci1=ci[1],ci2=ci[2])
		
			write.table(data_row, file = outFile, sep = ",", append = TRUE, quote = FALSE,  col.names = FALSE, row.names = FALSE)
		}
	}
	
	
}
