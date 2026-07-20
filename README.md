# Benchmark for Kronecker Product with Custom User Selector

Данный репозиторий содержит бенчмарки для демонстрации работы селектора в произведении Кронекера в библиотеке SuiteSparse:GraphBLAS

## Структура
- `GraphBLAS/` --- сабмодуль с модифицированной версией библиотеки (ветка `feat/kroner-user-selector`)
- `src/` --- исходный код бенчмарков (JIT/без JIT версии, с заданным пользователем/встроенным селектором и без него)
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
Селектор --- это стандартная для селектора функция типа `GrB_IndexUnaryOp` со следующей сигнатурой:
```bash
void selector(bool *result, const void *value, GrB_Index i, GrB_Index j, const void *y);
```
Параменты:
- `result` - указатель на bool: true - сохранение элемента, false - отбрасывание;
- `value` - указатель на значение элемента матрицы;
- `i, j` - индексы строки и столбца элемента в результирующей матрице;
- `y` - опциональный скалярный параметр.

Если селектор равен `NULL`, поведение по умолчанию меняется: фильтрация не применяется, и все вычисленные значения сохраняются в результат.