from pathlib import *
fna = PurePath(
    '/home/cszhang/workSpace/classification/refSeq/pipline/for_kraken_db/library.fna')
out = fna.parent.joinpath('{file_name}.removeNX{suffix}'.format(
    file_name=fna.stem, suffix=fna.suffix))
out_steam = open(out, 'w')
with open(fna, 'r') as f:
    line = f.readline()
    while line:
        if not line.startswith('>'):
            line = line.replace('n', '').replace(
                'N', '').replace('x', '').replace('X', '')
        out_steam.write(line)
        line = f.readline()
out_steam.flush()
out_steam.close()
