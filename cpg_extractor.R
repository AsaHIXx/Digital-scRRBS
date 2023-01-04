## special for methylkit CpG extractor / input: sorted bam file
library('getopt')
library('methylKit')
command=matrix(c(
    'help', 'h', 0,'loical', '显示此帮助信息',
    'input', 'i', 1, 'character', '输入文件',
    'path', 'o', 1, 'character', '输出位置'),
    byrow=T, ncol=5
)
file_name <- read.table(command$input, sep = '\t', names=c('file'))
file <- as.list(file_name$file)
file_num <- c()
for (i in as.character(file)){
    file_num <- append(file_num, strsplit(i, '_')[[1]][1])
}
file_num <- as.list(file_num)
my.methRaw=processBismarkAln(location=file,
                         sample.id=file_num, assembly="hg19",
                        treatment=as.double(rep(1, times=length(file))),mincov=0,
                         read.context=c("CpG"), save.folder=command$path)