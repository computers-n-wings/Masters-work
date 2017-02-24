void Print_Matrix(double matrix[], int n)
{
	for (int i = 0; i<n; ++i)
	{
		for (int j = 0; j<n; ++j)
		{
			cout << setw(10) << matrix[i*n + j] << setw(10);
		} 
		cout << endl;
	}
}

void Build_K_elemental(double matrix[], int n, double a, double b, double c)
{
	for (unsigned int i = 0; i<n; ++i)
	{
		for (unsigned int j = 0; j<n; ++j)
		{
			if (i <= j)
			{
				if (i == 0)
				{
					Ke[i*n+0] = a;
					Ke[i*n+3] = -a;
					if (j!=0 && j!=3)
					{
						Ke[i*n+j] = 0;
					}
				}
				if (i == 1)
				{
					Ke[i*n+1] = b;
					Ke[i*n+2] = c;
					Ke[i*n+5] = c;
					Ke[i*n+4] = -b;
					Ke[i*n+3] = 0;
				}
				if (i == 2)
				{
					Ke[i*n+2] = 4*a;
					Ke[i*n+4] = -c;
					Ke[i*n+5] = 2*a;
					Ke[i*n+3] = 0;
				}
				if (i == 3)
				{
					Ke[i*n+3] = a;
					if (j!=3)
					{
						Ke[i*n+j] = 0;
					}
				}
				if (i == 4)
				{
					Ke[i*n+4] = b;
					Ke[i*n+5] = -c;
				}
				if (i == 5)
				{
					Ke[i*n+5] = 4*a;
				}
			}
			else
			{
				Ke[i*n+j] = Ke[j*n+i];
			}
		}
	}
}