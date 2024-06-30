#include <vector>
#include <math.h>

// *** QMC�p�֐��� ***
std::vector<int> prime_numbers;

void init_prime_numbers() {
	// �G���g�X�e�l�X���
	const int N = 100000000; // 10000000�܂ł̑f�����e�[�u���Ɋi�[����i664579�j
	const int sqrtN = sqrt((double)N);
	std::vector<char> table;
	table.resize(N + 1, 0);

	for (int i = 2; i*i <= N; i++) {
		if (table[i] == 0) { // i���f���ł�����
			// �ӂ邤
			for (int j = i + i; j <= N; j += i)
				table[j] = 1;
		}
	}

	// �񋓂���
	for (int i = 2; i <= N; i++) {
		if (table[i] == 0) {
			prime_numbers.push_back(i);
		}
	}
	fprintf(stderr, "prime numbers:%d\n", prime_numbers.size());
}

