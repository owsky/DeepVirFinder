$start = Get-Date -UFormat "%s.%N"
$lengths = @(150, 300, 500, 1000)

foreach ($l in $lengths) {
  python encode.py -i ./train_data/tr/host_tr.fa -l $l -p host
  python encode.py -i ./train_data/tr/virus_tr.fa -l $l -p virus
  python encode.py -i ./train_data/val/host_val.fa -l $l -p host
  python encode.py -i ./train_data/val/virus_val.fa -l $l -p virus
}

$end = Get-Date -UFormat "%s.%N"
$runtime = [decimal]($end - $start)
Write-Host "Running time for encoding is $runtime seconds"

foreach ($l in $lengths) {
  python training.py -l $l -i ./train_data/tr/encode -j ./train_data/val/encode -o ./models -f 10 -n 500 -d 500 -e 10
  $end = Get-Date -UFormat "%s.%N"
  $runtime = [decimal]($end - $start)
  Write-Host "Running time for training with length $l is $runtime seconds"
}

$end = Get-Date -UFormat "%s.%N"
$runtime = [decimal]($end - $start)
Write-Host "Global running time is $runtime seconds"
