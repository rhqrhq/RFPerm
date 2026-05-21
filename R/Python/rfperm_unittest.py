#unittes:


def _validate_loss_output(loss_values, expected_n, name = 'loss'):
    loss_values = np.asarray(loss_values, dtype = float).reshape(-1)
    if loss_values.shape[0] != expected_n:
        raise ValueError(
            f"{name} must return one loss value per observation. "
            f"Expected length {expected_n}, got {loss_values.shape[0]}."
        	)
    if not np.all(np.isfinite(loss_values)):
        raise ValueError(f"{name} contains NaN or infinite values.")
    return loss_values


def _validate_rfperm_inputs(df_exist, df_new, loss, alpha, B, n_perm):
    if not isinstance(df_exist, pd.DataFrame):
        raise TypeError('df_exist must be a pandas DataFrame')
    if not isinstance(df_new, pd.DataFrame):
        raise TypeError('df_new must be a pandas DataFrame')
    if df_exist.shape[1] < 2 or df_new.shape[1] < 2:
        raise ValueError('df_exist and df_new must contain at least one feature column and one response column.')
    if df_exist.shape[1] != df_new.shape[1]:
        raise ValueError(
            f"df_exist and df_new must have the same number of columns. "
            f"God {df_exist.shape[1]} and {df_new.shape[1]}."
        )
    if list(df_exist.columns) != list(df_new.columns):
        raise ValueError('df_exist and df_new must have the same column names in the same order.')
    if df_exist.isna().any().any() or df_new.isna().any().any():
        raise ValueError('df_exist and df_new contain missing values. Please remove or impute them first.')
    if not callable(loss):
        raise TypeError("loss must be callable.")
    if not (0 < alpha < 1):
        raise ValueError('alpha must be between 0 and 1.')
    if not instance(B, int) or B <= 0:
        raise ValueError('B must be a positive integer.')
    if not instance(n_perm, int) or n_perm <= 0:
        raise ValueError('n_perm must be a positive integer.')
    return 











