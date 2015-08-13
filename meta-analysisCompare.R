# Meta-analysis comparison

library(MAMA)
data(ColonData)

pval<- metaMA(ColonData, "MSIstatus", which ="pval")

es <- metaMA(ColonData, "MSIstatus", which = "ES")

SOGL.res <- performSOGL(ColonData, varname = "MSIstatus",test = "FCH", B = 100, which=c("score", "empirical"),min.weight = 1e-05, two.sided = TRUE )

rp <- RankProduct(ColonData, "MSIstatus", num.perm = 50,gene.names = rownames(GEDM(ColonData)[[1]]) )

z.stat<-posterior.mean(ColonData, "MSIstatus", nsamp=5)

vm <- VennMapper(ColonData, varname="MSIstatus", cutoff=1)

map <- MAP.Matches(ColonData, "MSIstatus", t.cutoff = "95.00%",nperm = 100, sig.col = "p.col.strong")

metra<-METRADISC(ColonData, "MSIstatus", nperm = 1000)

lists<-join.DEG(pval, es, SOGL.res, rp, z.stat,map, cutoff=0.05)
names(lists)<-c("P-value", "Effect size","Ordered List","RankProduct", "Z-statistics","MAP match")
summary(lists)

MAT<-make.matrix(lists)
MAT[1:5,1:5]

n.met<-apply(MAT,1,sum)
hist(n.met, main="", xlab="Number of methods",ylab="Number of genes", xlim=c(1,8))

dim(MAT[n.met>5,])
# [1] 48  7
n.gen<-apply(MAT,2,sum)


p <- barplot(n.gen, cex.names=0.8, las=2, ylab="Number of genes", col = c("lightblue","mistyrose", "lightcyan", "lavender", "cornsilk","blue"))

pdf("MetaComparison.pdf", width=6, height=6)
print(p)
dev.off()