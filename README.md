# nbody — физическое ядро

2D-симуляция N тел на Go: прямой N², Барнс—Хат, leapfrog.

    go test ./tests/ -v          # тесты
    go run ./cmd/sim -scenario dense_core -n 2000 -steps 600

Сценарии: uniform, dense_core, collision. Решатели: `-solver direct|bh`, `-theta 0.5`, `-workers N`.
