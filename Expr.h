#ifndef __Expr__
#define __Expr__

#include <vector>

#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <string>


/* error codes */
enum EvalErrors
{
	EXPR_GOOD,		/* expression totally good */
	EXPR_SOSO,		/* expression partially good */
	EXPR_BAD		/* expression totally bad */
};


double expr_eval(const char* str, std::vector<double*>& value);
double expr_eval(const char* str, std::vector<double>& value);
double	expr_eval(const char* s);

int set_Variable(std::string name, double value);

#define M_PI       3.14159265358979323846

class Variable
{
public:
	std::string name;
	double value;
};
extern std::vector<Variable> VariableList;

class Expr
{
	char *s0, *s;
	int expr_error;

	inline double DegsToRads(double degrees)	///< Degrees to radians conversion...
	{
		return(degrees * (M_PI / 180.0));
	}
	inline double RadsToDegs(double rads)		///< Radians to degrees conversion...
	{
		return(rads * (180.0 / M_PI));
	}

	inline void space()
	{
		for (; isspace(*s); s++);
	}


	inline double term()
	{
		double x, y;

		for (x = factor();;)
		{
			//space();
			switch (*s)
			{
			case '*':
				s++;
				x *= factor(); break;
			case '/': s++; x /= factor(); break;
			case '%': s++; y = factor(); x = x - floor(x / y) * y; break;
			default: return x;
			}
		}
	}

	inline double factor()
	{
		double x;

		for (x = signednumber();;)
		{
			//space();
			switch (*s)
			{
			case '^':
				s++;
				return pow(x, factor());  	/* right-associative */
			default: return x;
			}
		}
	}

	inline double signednumber()
	{
		space();
		switch (*s)
		{
		case '-':
			s++;
			return -signednumber();
		case '+': s++; return signednumber();
		default: return number();
		}
	}

	inline double number()
	{
		const char *func;
		int n;
		double x, y;

		space();
		if (isdigit(*s) || *s == '.') return posconst();
		if (*s == '(') return paren();

		if (isalpha(*s))
		{
			func = s;
			for (s++; isalpha(*s) || isdigit(*s); s++);
			n = s - func;  	/* length of funcname */

			if (eq(n, func, "pi"))	return M_PI;
			if (eq(n, func, "e"))	return exp(1.);

			if (eq(n, func, "sqrt"))	return sqrt(paren());
			if (eq(n, func, "exp"))	return exp(paren());
			if (eq(n, func, "log"))	return log(paren());
			if (eq(n, func, "pow"))
			{
				paren2(&x, &y);
				return pow(x, y);
			}

			if (eq(n, func, "sin"))	return sin(paren());
			if (eq(n, func, "cos"))	return cos(paren());
			if (eq(n, func, "tan"))	return tan(paren());
			if (eq(n, func, "asin"))	return asin(paren());
			if (eq(n, func, "acos"))	return acos(paren());
			if (eq(n, func, "atan"))	return atan(paren());
			if (eq(n, func, "atan2"))
			{
				paren2(&x, &y);
				return atan2(x, y);
			}

			if (eq(n, func, "sind"))	return sin(DegsToRads(paren()));
			if (eq(n, func, "cosd"))	return cos(DegsToRads(paren()));
			if (eq(n, func, "tand"))	return tan(DegsToRads(paren()));
			if (eq(n, func, "dasin"))	return RadsToDegs(asin(paren()));
			if (eq(n, func, "dacos"))	return RadsToDegs(acos(paren()));
			if (eq(n, func, "datan"))	return RadsToDegs(atan(paren()));
			if (eq(n, func, "datan2"))
			{
				paren2(&x, &y);
				return RadsToDegs(atan2(x, y));
			}

			if (eq(n, func, "floor"))	return floor(paren());
			if (eq(n, func, "ceil"))	return ceil(paren());

			for (int i = 0; i < VariableList.size(); i++)
			{
				if (eq(n, func, VariableList[i].name.c_str()))	return VariableList[i].value;
			}
			error(func, n, "bad numerical expression");
			return 0.;
		}

		error(s, 1, "syntax error");
		return 0.;
	}

	/* paren: '(' expr ')' */

	inline double paren()
	{
		double x;

		space();
		if (*s != '(') error(s, 1, "expected '('");
		s++;
		x = expr();
		space();
		if (*s != ')') error(s, 1, "expected ')'");
		s++;
		return x;
	}

	/* paren2: '(' expr ',' expr ')' */

	inline void paren2(double *x, double *y)
	{
		space();
		if (*s != '(') error(s, 1, "expected '('");
		s++;
		*x = expr();
		space();
		if (*s != ',') error(s, 1, "expected ','");
		s++;
		*y = expr();
		space();
		if (*s != ')') error(s, 1, "expected ')'");
		s++;
	}

	/*
	* posconst: given a std::string& beginning at s, return floating point value.
	* like atof but it uses and modifies the global ptr s
	*/

	inline double posconst()
	{
		int base, exp, pos, d;
		double x, y;

		space();
		if (*s == '0')
		{		/* change base: 10 = 012 = 0xa = 0b2:1010 */
			s++;
			switch (*s)
			{
			case 'b':
				s++;
				for (base = 0; isdigit(*s); s++)
					base = base * 10 + *s - '0';  	/* base is in base 10! */
				if (*s != ':') error(s, 1, "expecting ':'");
				s++;
				break;
			case 'x':
				s++;
				base = 16;
				break;
			case 't':
				s++;
				base = 10;
				break;
			case '.':
				base = 10;
				break;  		/* a float, e.g.: 0.123 */
			default:
				base = 8;
				break;
			}
		}
		else base = 10;

		x = 0.;
		for (; d = digit(*s), d >= 0 && d < base; s++)
			x = x * base + d;
		if (*s == '.')
		{
			s++;
			for (y = 1.; d = digit(*s), d >= 0 && d < base; s++)
			{
				x = x * base + d;  		/* fraction is in variable base */
				y *= base;
			}
			x /= y;
		}
		if (*s == 'e' || *s == 'E')
		{
			s++;
			if (*s == '-')
			{
				s++;
				pos = 0;
			}
			else if (*s == '+')
			{
				s++;
				pos = 1;
			}
			else pos = 1;
			for (exp = 0; isdigit(*s); s++)
				exp = exp * 10 + *s - '0';  	/* exponent is in base 10 */
			y = expt(base, exp);
			if (pos) x *= y;
			else x /= y;
		}
		return x;
	}

	inline int digit(int c)
	{
		return isdigit(c) ? c - '0' :
			c >= 'a' && c <= 'z' ? c - 'a' + 10 : c >= 'A' && c <= 'Z' ? c - 'A' + 10 : -1;
	}

	/* expt: a^n for n>=0 */

	inline double expt(int a, int n)
	{
		double t, x;

		if (n < 0)
		{
			fprintf(stderr, "expt: can't do negative exponents\n");
			return 1.;
		}
		if (n == 0) return 1.;
		for (t = a, x = 1.; n > 0; n >>= 1)
		{
			if (n & 1) x *= t;
			t *= t;
		}
		return x;
	}

	/* eq: test equality of std::string& a, of length n, with null-terminated std::string& b */

	inline int eq(int n, const char *a, const char *b)
	{
		char c, *nca = (char *)a;
		int ret;

		c = a[n];
		nca[n] = 0;
		ret = (strcmp(a, b) == 0);
		nca[n] = c;
		return ret;
	}


	/* prints: print std::string& s of length n */

	inline void prints(int n, const char *s)
	{
		char c, *ncs = (char *)s;

		c = s[n];
		ncs[n] = 0;
		printf("%s", s);
		ncs[n] = c;
	}

public:
	inline Expr()
	{
		s0 = s = 0;
		expr_error = 0;
	}
	inline Expr(char* exper)
	{
		if (strlen(exper) >= 2 && (exper[strlen(exper) - 2] == ' ' || exper[strlen(exper) - 2] == '\t'))
		{
			printf("Extra blank at the end [%s]\n", exper);
			fprintf(stderr, "Extra blank at the end [%s]\n", exper);
			fflush(stdout);
			fflush(stderr);
			throw "error expression.";
		}
		s0 = s = exper;
		expr_error = 0;
	}
	inline Expr(const char* exper)
	{
		if (strlen(exper) >= 2 && (exper[strlen(exper) - 2] == ' ' || exper[strlen(exper) - 2] == '\t'))
		{
			printf("Extra blank at the end [%s]\n", exper);
			fprintf(stderr, "Extra blank at the end [%s]\n", exper);
			fflush(stdout);
			fflush(stderr);
			throw "error expression.";
		}
		s0 = s = (char*)exper;
		expr_error = 0;
	}

	inline double expr()
	{
		double x;

		space();
		for (x = term();;)
		{
			//space();
			switch (*s)
			{
			case '+':
				s++;
				x += term(); break;
			case '-': s++; x -= term(); break;
			default: return x;
			}
		}
	}
	inline char* cur()
	{
		return s;
	}
	inline char* next()
	{
		s++;
		return s;
	}
	inline void error(const char *s, int len, const char *err)
	{
		//    if (*s == 0) s[len] = 0;  	/* just in case */
		printf("expr: %s: ", err);
		prints(s - s0, s0);
		printf("[");
		prints(len, s);
		printf("]");
		prints(s + strlen(s) - s0 - len, s + len);
		printf("\n");
		if (expr_error != EXPR_BAD)
			expr_error = s == s0 ? EXPR_BAD : EXPR_SOSO;
	}

	inline void error_chk()
	{
		expr_error = s == s0 ? EXPR_BAD : EXPR_SOSO;
	}
	inline int error()
	{
		return expr_error;
	}
};

#endif
