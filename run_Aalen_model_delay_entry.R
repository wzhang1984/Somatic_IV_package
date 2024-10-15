library(timereg)

train = read.csv('tmp/train.csv', row.names=1, check.names=F)
metadata = read.csv('tmp/metadata.csv', row.names=1, check.names=F)

res = list()
i = 0
for (gene in colnames(train)) {
    df = cbind(train[gene], metadata[rownames(train), ])
    df = df[complete.cases(df), ]
    colnames(df) = c('gene', colnames(metadata))
    x = 'const(gene)'
    for (confounder in tail(colnames(df), n=-4)) {
        x = paste0(x, '+const(', confounder, ')')
    }
    out = aalen(as.formula(paste('Surv(time=time_delay_entry, time2=time_to_event, event=event)', '~', x)), df, n.sim=100)
    out = coef.aalen(out)
    rownames(out) = c(gene)
    i = i+1
    res[[i]] = out
}

res = do.call(rbind, res)
write.table(res, 'tmp/results.txt', sep='\t', quote=F, col.names=NA)
