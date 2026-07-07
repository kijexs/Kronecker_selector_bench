# Benchmark for Kronecker Product with Custom User Selector

Данный репозиторий содержит бенчмарки для демонстрации работы селектора в произведении Кронекера в библиотеке SuiteSparse:GraphBLAS

## Структура
- `GraphBLAS/` --- сабмодуль с модифицированной версией библиотеки (ветка `feat/kroner-user-selector`)
- `src/` --- исходный код бенчмарков (JIT/без JIT версии, с заданным селектором и без)
- `kron_mat/` --- входные матрицы и метаданные
- `bench_res/` --- директория для сохраняемых результатов
- `e_meta.txt` --- количество запусков

## Подготовка данных
Входные данные (матрицы и метаданные) доступны по ссылке: https://drive.google.com/drive/folders/1fwZii07506EdDsPOFHTnGAnWYzvmJlG2

## Инструкция по сборке и запуску

### Клонирование репозитория
```bash
git clone --recurse-submodules https://github.com/kijexs/Kronecker_selector_bench.git
```

### Сборка графбласа
```bash
cd GraphBLAS

cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE="Debug" \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_C_COMPILER_LAUNCHER="ccache" \
    -DCMAKE_CXX_COMPILER_LAUNCHER="ccache" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cd build && ninja && cd ../..
```

### Сборка бенчмарков
```bash
gcc src/bench_jit_selector.c \
    -I./GraphBLAS/Include \
    -L./GraphBLAS/build \
    -lgraphblas -lm \
    -Wl,-rpath ./GraphBLAS/build \
    -o GraphBLAS/build/kron_bench
```
### Запуск
Пример для avrora_3
```bash
./GraphBLAS/build/kron_bench   kron_mat/rsa_avrora.csv   kron_mat/graph_avrora_3.csv   kron_mat/types_avrora.txt   e_meta.txt   bench_res/test_output.csv
```
В bench_res/test_output.csv находятся результаты замеров времени

Возможно перед этим необходимо будет выполнить:
```bash
export LD_LIBRARY_PATH=./GraphBLAS/build:$LD_LIBRARY_PATH
```

## Селектор
Селектор --- это функция типа `GxB_unary_function` со следующей сигнатурой:
```bash
typedef void (*GxB_unary_function)(void *, const void *)
```
Селектор вызывается как `sel(&result, &value)`, где result --- это bool, которому должно быть присвоено значение true, если значение должно быть сохранено, или false в противном случае

Если значение sel равно NULL, используется поведение по умолчанию: значение сохраняется, если какой-либо из его байтов отличен от нуля