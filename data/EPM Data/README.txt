***************************
DATA DESCRIPTION
***************************

1) my_epm1.csv

	This is the raw dataset as obtained from the github link 
	mentioned in Claudia's book.
---------------------------------------------------
2) my_epm1_processed.csv
	
	A processed dataframe from my_epm1.csv, that has the households as a single 
	column (instead of a pivot_wider format in which my_epm1.csv came in), 
	with their respective clusters they belong into.
	Also, training and testing division is indicated by a boolean column.
-------------------------------------------------------
3) quants-cluster<number>.csv

	A dataframe containing the estimated quantiles for cluster denoted by <number>.
	The columns are 
		Index, Household1, Household2, (etc.) and Quantile 
	The quantile goes from 0.5 to 0.99.


