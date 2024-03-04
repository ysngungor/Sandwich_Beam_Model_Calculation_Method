# Code Scaffolding
+ In this study, I will focus on the code. You can find theoretical information on the subject in my master's thesis titled "An approach for the analysis of steel plate shear wall systems under earthquake loads", published in 2022.
[Thesis Center](https://tez.yok.gov.tr/UlusalTezMerkezi/giris.jsp)

+ After entering various data into the program, I performed simple mathematical operations with theoretical formulations.

+ I collected the data I received with the array matrix command in a table by creating a dataframe.

+ If the index exists in this table, the directly corresponding values are used.

```
y = False
for i in index_k:
    if i == k:
        df = df.loc[k]
        y = True
```

+ If it is not in the table, I interpolated a value between two rows in the table to find the corresponding values in other columns.

```
if y != True:
    df.loc[k] = np.nan  # Adding unknown row based on k value
    df = df.reindex(sorted(df.index), axis=0)  # Reordering the list with k value
    df = df.interpolate(method="polynomial", order=1).round(4)  # 1st degree polynomial interpolation
    df = df.loc[k]
```

### For Example;

1. Let's have a table like below.

|index|S1|S2|S3|
|---|---|---|---|
|1.2|2|5|11|
|1.8|4|7|5|
|2.1|6|8|3|
</br>

2. Even if they do not appear in this table, if we want to find the values ​​corresponding to the index value 1.6, we add a new row with the values nan.

|index|S1|S2|S3|
|---|---|---|---|
|1.2|2|5|11|
|1.8|4|7|5|
|2.1|6|8|3|
|***1.6***|***none***|***none***|***none***|
</br>

3. We reorder the table according to this index value.

|index|S1|S2|S3|
|---|---|---|---|
|1.2|2|5|11|
|***1.6***|***none***|***none***|***none***|
|1.8|4|7|5|
|2.1|6|8|3|
</br>

4. Then, by performing polynomial interpolation, we find the values of S1, S2 and S3 corresponding to the value 1.6.
Since the value 1.6 is closer to 1.8, we expect the other column values to be closer to the values corresponding to 1.8.

|index|S1|S2|S3|
|---|---|---|---|
|1.2|2|5|11|
|***1.6***|***3.33***|***6.33***|***7***|
|1.8|4|7|5|
|2.1|6|8|3|
</br>

+ After then, mathematical operations were performed according to the data obtained from the general data table and thus the final result was reached.


