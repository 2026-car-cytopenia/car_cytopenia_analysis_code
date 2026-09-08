# Shared settings helper

`settings.R` supplies two functions to the original clinical/single-cell scripts:
`study_path()` resolves input/output paths against a private configuration file;
`study_setting()` reads sample selections, crosswalks and display choices.
Each analysis directory has its own `settings.example.R`. This helper is loaded
by those scripts; it is not a separate analysis or a pipeline runner.
