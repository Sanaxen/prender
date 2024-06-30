#include <vector>
#include <math.h>

// *** QMC用関数等 ***
std::vector<int> prime_numbers;

void init_prime_numbers() {
	// エラトステネスの篩
	const int N = 100000000; // 10000000までの素数をテーブルに格納する（664579個）
	const int sqrtN = sqrt((double)N);
	std::vector<char> table;
	table.resize(N + 1, 0);

	for (int i = 2; i*i <= N; i++) {
		if (table[i] == 0) { // iが素数であった
			// ふるう
			for (int j = i + i; j <= N; j += i)
				table[j] = 1;
		}
	}

	// 列挙する
	for (int i = 2; i <= N; i++) {
		if (table[i] == 0) {
			prime_numbers.push_back(i);
		}
	}
	fprintf(stderr, "prime numbers:%d\n", prime_numbers.size());
}

