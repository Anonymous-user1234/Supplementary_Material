
class Enumerator {

private:

	wind_t len = 0; // the number of indexes
	wind_t range = 0; // the max valid position is range - 1;

	std::vector<wind_t> sel_idxs;
	wind_t curr_pos = 0;

	bool first = true;
	wind_t ctr = 0;

public:

	Enumerator(wind_t _len, wind_t _range) {
		len = _len;
		range = _range;

		assert(len > 0);
		sel_idxs.clear();
		sel_idxs.resize(len);
		sel_idxs[0] = 0;
		
		curr_pos = 0;
		assert(range >= 1);

		first = true;
		ctr = 0;
	}

	std::vector<wind_t> const & current() {
		return sel_idxs;
	}

	void print() {
		std::string st = "indexes = [";
		for (auto & idx : sel_idxs) {
			if (st != "indexes = [")
				st += " ";
			st += std::to_string(idx);
		}
		st += "], idx = " + std::to_string(ctr);

		std::cout << std::endl;
		PRINTF_STAMP("%s\n", st.c_str());
	}

	bool step() { //ok4.6

		if (!first) {
			assert(curr_pos == len - 1);
			++sel_idxs[curr_pos];
		}

		while (true) {
			if (sel_idxs[curr_pos] >= range) {
				--curr_pos;
				if (curr_pos < 0) break;
				++sel_idxs[curr_pos];
				continue;
			}

			if (curr_pos < len - 1) {
				// move forward
				sel_idxs[curr_pos + 1] = sel_idxs[curr_pos] + 1;
				++curr_pos;
				continue;
			}

			assert(curr_pos == len - 1);

			if (first) ctr = 0;
			else ++ctr;

			first = false;
			
			return true;
			
			// ++sel_idxs[curr_pos];
		}

		return false;
	}


};
