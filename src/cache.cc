#include "cache.h"
#include "set.h"
#include "block.h"
#include <algorithm>
#include "ooo_cpu.h"

#define DEBUG_STLB 0
#define DEBUG 0


void CACHE::print_fctb(){
	for(int i=0; i<FCTB_SIZE; i++)
		cout << hex << fctb[i][0] << ", " <<  dec << fctb[i][1] << ", " << fctb[i][2] << ", " << fctb[i][3] << endl;
}

int CACHE::search_fctb(uint64_t current_vpn){
	int fctb_found_pos = -10, acc;
	for(int f=0; f<FCTB_SIZE; f++){
		if(fctb[f][0] == current_vpn)
			return f;

		for(int t=1; t<=fctb[f][1]; t++){
			if((fctb[f][0]-t) == current_vpn)
				return f;
		}

		acc = 1;
		for(int t=fctb[f][1]+1; t<8; t++){
			if((fctb[f][0]+acc) == current_vpn)
				return f;
			acc++;
		}
	}
	return -10;
}


int CACHE::fctb_replacement_policy(){
	uint64_t lru_min = fctb[0][2];
	int lru_victim = 0;
	for(int f=1; f<FCTB_SIZE; f++){
		if(fctb[f][2] < lru_min){
			lru_min = fctb[f][2];
			lru_victim = f;
		}
	}
	return lru_victim;
}

void CACHE::refresh_fctb(uint64_t current_cycle){
	for(int f=0; f<FCTB_SIZE; f++){
		//if((current_cycle - fctb[f][2]) > PAGE_TABLE_LATENCY)
		if(current_cycle > 1*(fctb[f][2] + fctb[f][3])){
			fctb[f][0] = 0;
			fctb[f][1] = 0;
			fctb[f][2] = 0;
			fctb[f][3] = 0;
		}
	}
}

int * CACHE::sorted_free_distances(){
	/*
	   uint64_t sorted_table[14];

	   for(int i=0; i<14; i++)
	   sorted_table[i] = free_distance_table[i];

	   int n = sizeof(sorted_table)/sizeof(sorted_table[0]);

	   sort(sorted_table, sorted_table+n);

	//for(int i=0; i<14; i++){
	//    cout << sorted_table[i] << ", ";
	//}
	//cout << endl;
	*/

	/* pair sort */
	int n = sizeof(free_distance_table)/sizeof(free_distance_table[0]);
	pair<uint64_t, int> * pairf;
	pairf = (pair<uint64_t, int>*) malloc(14*sizeof( pair<uint64_t, int>));

	for(int i=0; i<14; i++){
		pairf[i].second = i - 6;
		if(i <= 6)
			pairf[i].second--;
		pairf[i].first = free_distance_table[i];
	}

	sort(pairf, pairf+n);

	/*
	   for(int i=0; i<14; i++){
	   cout << pairf[i].first << ", ";
	   }
	   cout << endl;
	   for(int i=0; i<14; i++){
	   cout << pairf[i].second << ", ";
	   }
	   */

	int * indexes;
	indexes = (int *) malloc(14*sizeof(int));
	for(int i=0; i<14; i++){
		indexes[i] = pairf[i].second;
	}

	free(pairf);

	return indexes;

	//return pairf;
}



void CACHE::print_pml4(){
	for(int s=0; s<PML4_SET; s++){
		for(int w=0; w<PML4_WAY; w++)
			cout << pml4[s][w] << ", " << pml4_lru[s][w];
		cout << endl;
	}
}

void CACHE::print_pdp(){
	for(int s=0; s<PDP_SET; s++){
		for(int w=0; w<PDP_WAY; w++)
			cout << pdp[s][w] << ", " << pdp_lru[s][w];
		cout << endl;
	}
}

void CACHE::print_pd(){
	for(int s=0; s<PD_SET; s++){
		for(int w=0; w<PD_WAY; w++)
			cout << pd[s][w] << ", " << pd_lru[s][w] << " ;; ";
		cout << endl;
	}
}

int CACHE::search_pml4(uint64_t address){
	for(int s=0; s<PML4_SET; s++){
		for(int w=0; w<PML4_WAY; w++){
			if(pml4[s][w] == address)
				return 1;
		}
	}
	return 0;
}

int CACHE::search_pdp(uint64_t address){
	for(int s=0; s<PDP_SET; s++){
		for(int w=0; w<PDP_WAY; w++){
			if(pdp[s][w] == address)
				return 1;
		}
	}
	return 0;
}

int CACHE::search_pd(uint64_t address){
	for(int s=0; s<PD_SET; s++){
		for(int w=0; w<PD_WAY; w++){
			if(pd[s][w] == address)
				return 1;
		}
	}
	return 0;
}

void CACHE::lru_pml4(uint64_t timer, uint64_t address){
	for(int s=0; s<PML4_SET; s++){
		for(int w=0; w<PML4_WAY; w++){
			if(pml4[s][w] == address){
				pml4_lru[s][w] = timer;
				return;
			}
		}
	}

	uint64_t min = pml4_lru[0][0];
	int victim[2] = {0,0};

	for(int s=0; s<PML4_SET; s++){
		for(int w=0; w<PML4_WAY; w++){
			if(pml4_lru[s][w] < min){
				min = pml4_lru[s][w];
				victim[0] = s;
				victim[1] = w;
			}
		}
	}

	pml4[victim[0]][victim[1]] = address;
	pml4_lru[victim[0]][victim[1]] = timer;
	return;
}
void CACHE::lru_pdp(uint64_t timer, uint64_t address){
	for(int s=0; s<PDP_SET; s++){
		for(int w=0; w<PDP_WAY; w++){
			if(pdp[s][w] == address){
				pdp_lru[s][w] = timer;
				return;
			}
		}
	}

	uint64_t min = pdp_lru[0][0];
	int victim[2] = {0,0};

	for(int s=0; s<PDP_SET; s++){
		for(int w=0; w<PDP_WAY; w++){
			if(pdp_lru[s][w] < min){
				min = pdp_lru[s][w];
				victim[0] = s;
				victim[1] = w;
			}
		}
	}

	pdp[victim[0]][victim[1]]      = address;
	pdp_lru[victim[0]][victim[1]]  = timer;
	return;
}
void CACHE::lru_pd(uint64_t timer, uint64_t address){
	int bits = ceil(log2(PD_SET));
	int index = address & ((1<<bits)-1);
	//cout << "Index: " << index << endl;

	for(int w=0; w<PD_WAY; w++){
		if(pd[index][w] == address){
			pd_lru[index][w] = timer;
			return;
		}
	}

	uint64_t min = pd_lru[index][0];
	int victim[2] = {index,0};

	for(int w=0; w<PD_WAY; w++){
		if(pd_lru[index][w] < min){
			min = pd_lru[index][w];
			victim[0] = index;
			victim[1] = w;
		}
	}

	pd[victim[0]][victim[1]]      = address;
	pd_lru[victim[0]][victim[1]]  = timer;

	return;
}


int compare (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}


uint64_t l2pf_access = 0;

void CACHE::handle_fill()
{
	// handle fill
	uint32_t fill_cpu = (MSHR.next_fill_index == MSHR_SIZE) ? NUM_CPUS : MSHR.entry[MSHR.next_fill_index].cpu;
	if (fill_cpu == NUM_CPUS)
		return;

	if (MSHR.next_fill_cycle <= current_core_cycle[fill_cpu]) {

#ifdef SANITY_CHECK
		if (MSHR.next_fill_index >= MSHR.SIZE)
			assert(0);
#endif

		uint32_t mshr_index = MSHR.next_fill_index;

		// find victim
		uint32_t set = get_set(MSHR.entry[mshr_index].address), way;
		if (cache_type == IS_LLC) {
			way = llc_find_victim(fill_cpu, MSHR.entry[mshr_index].instr_id, set, block[set], MSHR.entry[mshr_index].ip, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].type);
		}
		else
			way = find_victim(fill_cpu, MSHR.entry[mshr_index].instr_id, set, block[set], MSHR.entry[mshr_index].ip, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].type);

#ifdef LLC_BYPASS
		if ((cache_type == IS_LLC) && (way == LLC_WAY)) { // this is a bypass that does not fill the LLC

			// update replacement policy
			if (cache_type == IS_LLC) {
				llc_update_replacement_state(fill_cpu, set, way, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].ip, 0, MSHR.entry[mshr_index].type, 0);

			}
			else
				update_replacement_state(fill_cpu, set, way, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].ip, 0, MSHR.entry[mshr_index].type, 0);

			// COLLECT STATS
			sim_miss[fill_cpu][MSHR.entry[mshr_index].type]++;
			sim_access[fill_cpu][MSHR.entry[mshr_index].type]++;

			// check fill level
			if (MSHR.entry[mshr_index].fill_level < fill_level) {

				if(fill_level == FILL_L2)
				{
					if(MSHR.entry[mshr_index].fill_l1i)
					{
						upper_level_icache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					}
					if(MSHR.entry[mshr_index].fill_l1d)
					{
						upper_level_dcache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					}
				}
				else
				{
					if (MSHR.entry[mshr_index].instruction)
						upper_level_icache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					if (MSHR.entry[mshr_index].is_data)
						upper_level_dcache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
				}
			}

			if(warmup_complete[fill_cpu] && (MSHR.entry[mshr_index].cycle_enqueued != 0))
			{
				uint64_t current_miss_latency = (current_core_cycle[fill_cpu] - MSHR.entry[mshr_index].cycle_enqueued);
				total_miss_latency += current_miss_latency;
			}

			MSHR.remove_queue(&MSHR.entry[mshr_index]);
			MSHR.num_returned--;

			update_fill_cycle();

			return; // return here, no need to process further in this function
		}
#endif

		uint8_t  do_fill = 1;

		// is this dirty?
		if (block[set][way].dirty) {

			// check if the lower level WQ has enough room to keep this writeback request
			if (lower_level) {
				if (lower_level->get_occupancy(2, block[set][way].address) == lower_level->get_size(2, block[set][way].address)) {

					// lower level WQ is full, cannot replace this victim
					do_fill = 0;
					lower_level->increment_WQ_FULL(block[set][way].address);
					STALL[MSHR.entry[mshr_index].type]++;

					DP ( if (warmup_complete[fill_cpu]) {
							cout << "[" << NAME << "] " << __func__ << "do_fill: " << +do_fill;
							cout << " lower level wq is full!" << " fill_addr: " << hex << MSHR.entry[mshr_index].address;
							cout << " victim_addr: " << block[set][way].tag << dec << endl; });
				}
				else {
					PACKET writeback_packet;

					writeback_packet.fill_level = fill_level << 1;
					writeback_packet.cpu = fill_cpu;
					writeback_packet.address = block[set][way].address;
					writeback_packet.full_addr = block[set][way].full_addr;
					writeback_packet.data = block[set][way].data;
					writeback_packet.instr_id = MSHR.entry[mshr_index].instr_id;
					writeback_packet.ip = 0; // writeback does not have ip
					writeback_packet.type = WRITEBACK;
					writeback_packet.event_cycle = current_core_cycle[fill_cpu];

					lower_level->add_wq(&writeback_packet);
				}
			}
#ifdef SANITY_CHECK
			else {
				// sanity check
				if (cache_type != IS_STLB)
					assert(0);
			}
#endif
		}

		if (do_fill){
			// update prefetcher
			if (cache_type == IS_L1I)
				l1i_prefetcher_cache_fill(fill_cpu, ((MSHR.entry[mshr_index].ip)>>LOG2_BLOCK_SIZE)<<LOG2_BLOCK_SIZE, set, way, (MSHR.entry[mshr_index].type == PREFETCH) ? 1 : 0, ((block[set][way].ip)>>LOG2_BLOCK_SIZE)<<LOG2_BLOCK_SIZE);
			if (cache_type == IS_L1D)
				l1d_prefetcher_cache_fill(MSHR.entry[mshr_index].full_addr, set, way, (MSHR.entry[mshr_index].type == PREFETCH) ? 1 : 0, block[set][way].address<<LOG2_BLOCK_SIZE,
						MSHR.entry[mshr_index].pf_metadata);
			if  (cache_type == IS_L2C)
				MSHR.entry[mshr_index].pf_metadata = l2c_prefetcher_cache_fill(MSHR.entry[mshr_index].address<<LOG2_BLOCK_SIZE, set, way, (MSHR.entry[mshr_index].type == PREFETCH) ? 1 : 0,
						block[set][way].address<<LOG2_BLOCK_SIZE, MSHR.entry[mshr_index].pf_metadata);
			if (cache_type == IS_LLC)
			{
				cpu = fill_cpu;
				MSHR.entry[mshr_index].pf_metadata = llc_prefetcher_cache_fill(MSHR.entry[mshr_index].address<<LOG2_BLOCK_SIZE, set, way, (MSHR.entry[mshr_index].type == PREFETCH) ? 1 : 0,
						block[set][way].address<<LOG2_BLOCK_SIZE, MSHR.entry[mshr_index].pf_metadata);
				cpu = 0;
			}

			// update replacement policy
			if (cache_type == IS_LLC) {
				llc_update_replacement_state(fill_cpu, set, way, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].ip, block[set][way].full_addr, MSHR.entry[mshr_index].type, 0);
			}
			else
				update_replacement_state(fill_cpu, set, way, MSHR.entry[mshr_index].full_addr, MSHR.entry[mshr_index].ip, block[set][way].full_addr, MSHR.entry[mshr_index].type, 0);

			// COLLECT STATS
			sim_miss[fill_cpu][MSHR.entry[mshr_index].type]++;
			sim_access[fill_cpu][MSHR.entry[mshr_index].type]++;

			fill_cache(set, way, &MSHR.entry[mshr_index]);

			// RFO marks cache line dirty
			if (cache_type == IS_L1D) {
				if (MSHR.entry[mshr_index].type == RFO)
					block[set][way].dirty = 1;
			}

			// check fill level
			if (MSHR.entry[mshr_index].fill_level < fill_level) {

				if(fill_level == FILL_L2)
				{
					if(MSHR.entry[mshr_index].fill_l1i)
					{
						upper_level_icache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					}
					if(MSHR.entry[mshr_index].fill_l1d)
					{
						upper_level_dcache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					}
				}
				else
				{
					if (MSHR.entry[mshr_index].instruction)
						upper_level_icache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
					if (MSHR.entry[mshr_index].is_data)
						upper_level_dcache[fill_cpu]->return_data(&MSHR.entry[mshr_index]);
				}
			}

			// update processed packets
			if (cache_type == IS_ITLB) { 
				MSHR.entry[mshr_index].instruction_pa = block[set][way].data;
				if (PROCESSED.occupancy < PROCESSED.SIZE)
					PROCESSED.add_queue(&MSHR.entry[mshr_index]);
			}
			else if (cache_type == IS_DTLB) {
				MSHR.entry[mshr_index].data_pa = block[set][way].data;
				if (PROCESSED.occupancy < PROCESSED.SIZE)
					PROCESSED.add_queue(&MSHR.entry[mshr_index]);
			}
			else if (cache_type == IS_L1I) {
				if (PROCESSED.occupancy < PROCESSED.SIZE)
					PROCESSED.add_queue(&MSHR.entry[mshr_index]);
			}
			//else if (cache_type == IS_L1D) {
			else if ((cache_type == IS_L1D) && (MSHR.entry[mshr_index].type != PREFETCH)) {
				if (PROCESSED.occupancy < PROCESSED.SIZE)
					PROCESSED.add_queue(&MSHR.entry[mshr_index]);
			}

			if(warmup_complete[fill_cpu] && (MSHR.entry[mshr_index].cycle_enqueued != 0))
			{
				uint64_t current_miss_latency = (current_core_cycle[fill_cpu] - MSHR.entry[mshr_index].cycle_enqueued);
				/*
				   if(cache_type == IS_L1D)
				   {
				   cout << current_core_cycle[fill_cpu] << " - " << MSHR.entry[mshr_index].cycle_enqueued << " = " << current_miss_latency << " MSHR index: " << mshr_index << endl;
				   }
				   */
				total_miss_latency += current_miss_latency;
			}

			MSHR.remove_queue(&MSHR.entry[mshr_index]);
			MSHR.num_returned--;

			update_fill_cycle();
		}
		}
	}

	void CACHE::handle_writeback()
	{
		// handle write
		uint32_t writeback_cpu = WQ.entry[WQ.head].cpu;
		if (writeback_cpu == NUM_CPUS)
			return;

		// handle the oldest entry
		if ((WQ.entry[WQ.head].event_cycle <= current_core_cycle[writeback_cpu]) && (WQ.occupancy > 0)) {
			int index = WQ.head;

			// access cache
			uint32_t set = get_set(WQ.entry[index].address);
			int way = check_hit(&WQ.entry[index]);

			if (way >= 0) { // writeback hit (or RFO hit for L1D)

				if (cache_type == IS_LLC) {
					llc_update_replacement_state(writeback_cpu, set, way, block[set][way].full_addr, WQ.entry[index].ip, 0, WQ.entry[index].type, 1);

				}
				else
					update_replacement_state(writeback_cpu, set, way, block[set][way].full_addr, WQ.entry[index].ip, 0, WQ.entry[index].type, 1);

				// COLLECT STATS
				sim_hit[writeback_cpu][WQ.entry[index].type]++;
				sim_access[writeback_cpu][WQ.entry[index].type]++;

				// mark dirty
				block[set][way].dirty = 1;

				if (cache_type == IS_ITLB)
					WQ.entry[index].instruction_pa = block[set][way].data;
				else if (cache_type == IS_DTLB)
					WQ.entry[index].data_pa = block[set][way].data;
				else if (cache_type == IS_STLB){
					WQ.entry[index].data = block[set][way].data;
					if(current_core_cycle[writeback_cpu] < block[set][way].stalls)
						add_stall_prefetch(block[set][way].stalls-current_core_cycle[writeback_cpu], writeback_cpu);
				}

				// check fill level
				if (WQ.entry[index].fill_level < fill_level) {

					if(fill_level == FILL_L2)
					{
						if(WQ.entry[index].fill_l1i)
						{
							upper_level_icache[writeback_cpu]->return_data(&WQ.entry[index]);
						}
						if(WQ.entry[index].fill_l1d)
						{
							upper_level_dcache[writeback_cpu]->return_data(&WQ.entry[index]);
						}
					}
					else
					{
						if (WQ.entry[index].instruction)
							upper_level_icache[writeback_cpu]->return_data(&WQ.entry[index]);
						if (WQ.entry[index].is_data)
							upper_level_dcache[writeback_cpu]->return_data(&WQ.entry[index]);
					}
				}

				HIT[WQ.entry[index].type]++;
				ACCESS[WQ.entry[index].type]++;

				// remove this entry from WQ
				WQ.remove_queue(&WQ.entry[index]);
			}
			else { // writeback miss (or RFO miss for L1D)

				DP ( if (warmup_complete[writeback_cpu]) {
						cout << "[" << NAME << "] " << __func__ << " type: " << +WQ.entry[index].type << " miss";
						cout << " instr_id: " << WQ.entry[index].instr_id << " address: " << hex << WQ.entry[index].address;
						cout << " full_addr: " << WQ.entry[index].full_addr << dec;
						cout << " cycle: " << WQ.entry[index].event_cycle << endl; });

				if (cache_type == IS_L1D) { // RFO miss

					// check mshr
					uint8_t miss_handled = 1;
					int mshr_index = check_mshr(&WQ.entry[index]);

					if(mshr_index == -2)
					{
						// this is a data/instruction collision in the MSHR, so we have to wait before we can allocate this miss
						miss_handled = 0;
					}
					else if ((mshr_index == -1) && (MSHR.occupancy < MSHR_SIZE)) { // this is a new miss

						if(cache_type == IS_LLC)
						{
							// check to make sure the DRAM RQ has room for this LLC RFO miss
							if (lower_level->get_occupancy(1, WQ.entry[index].address) == lower_level->get_size(1, WQ.entry[index].address))
							{
								miss_handled = 0;
							}
							else
							{
								add_mshr(&WQ.entry[index]);
								lower_level->add_rq(&WQ.entry[index]);
							}
						}
						else
						{
							// add it to mshr (RFO miss)
							add_mshr(&WQ.entry[index]);

							// add it to the next level's read queue
							//if (lower_level) // L1D always has a lower level cache
							lower_level->add_rq(&WQ.entry[index]);
						}
					}
					else {
						if ((mshr_index == -1) && (MSHR.occupancy == MSHR_SIZE)) { // not enough MSHR resource

							// cannot handle miss request until one of MSHRs is available
							miss_handled = 0;
							STALL[WQ.entry[index].type]++;
						}
						else if (mshr_index != -1) { // already in-flight miss

							// update fill_level
							if (WQ.entry[index].fill_level < MSHR.entry[mshr_index].fill_level)
								MSHR.entry[mshr_index].fill_level = WQ.entry[index].fill_level;

							if((WQ.entry[index].fill_l1i) && (MSHR.entry[mshr_index].fill_l1i != 1))
							{
								MSHR.entry[mshr_index].fill_l1i = 1;
							}
							if((WQ.entry[index].fill_l1d) && (MSHR.entry[mshr_index].fill_l1d != 1))
							{
								MSHR.entry[mshr_index].fill_l1d = 1;
							}

							// update request
							if (MSHR.entry[mshr_index].type == PREFETCH) {
								uint8_t  prior_returned = MSHR.entry[mshr_index].returned;
								uint64_t prior_event_cycle = MSHR.entry[mshr_index].event_cycle;
								MSHR.entry[mshr_index] = WQ.entry[index];

								// in case request is already returned, we should keep event_cycle and retunred variables
								MSHR.entry[mshr_index].returned = prior_returned;
								MSHR.entry[mshr_index].event_cycle = prior_event_cycle;
							}

							MSHR_MERGED[WQ.entry[index].type]++;

							DP ( if (warmup_complete[writeback_cpu]) {
									cout << "[" << NAME << "] " << __func__ << " mshr merged";
									cout << " instr_id: " << WQ.entry[index].instr_id << " prior_id: " << MSHR.entry[mshr_index].instr_id; 
									cout << " address: " << hex << WQ.entry[index].address;
									cout << " full_addr: " << WQ.entry[index].full_addr << dec;
									cout << " cycle: " << WQ.entry[index].event_cycle << endl; });
						}
						else { // WE SHOULD NOT REACH HERE
							cerr << "[" << NAME << "] MSHR errors" << endl;
							assert(0);
						}
					}

					if (miss_handled) {

						MISS[WQ.entry[index].type]++;
						ACCESS[WQ.entry[index].type]++;

						// remove this entry from WQ
						WQ.remove_queue(&WQ.entry[index]);
					}

				}
				else {
					// find victim
					uint32_t set = get_set(WQ.entry[index].address), way;
					if (cache_type == IS_LLC) {
						way = llc_find_victim(writeback_cpu, WQ.entry[index].instr_id, set, block[set], WQ.entry[index].ip, WQ.entry[index].full_addr, WQ.entry[index].type);
					}
					else
						way = find_victim(writeback_cpu, WQ.entry[index].instr_id, set, block[set], WQ.entry[index].ip, WQ.entry[index].full_addr, WQ.entry[index].type);

#ifdef LLC_BYPASS
					if ((cache_type == IS_LLC) && (way == LLC_WAY)) {
						cerr << "LLC bypassing for writebacks is not allowed!" << endl;
						assert(0);
					}
#endif

					uint8_t  do_fill = 1;

					// is this dirty?
					if (block[set][way].dirty) {

						// check if the lower level WQ has enough room to keep this writeback request
						if (lower_level) { 
							if (lower_level->get_occupancy(2, block[set][way].address) == lower_level->get_size(2, block[set][way].address)) {

								// lower level WQ is full, cannot replace this victim
								do_fill = 0;
								lower_level->increment_WQ_FULL(block[set][way].address);
								STALL[WQ.entry[index].type]++;

								DP ( if (warmup_complete[writeback_cpu]) {
										cout << "[" << NAME << "] " << __func__ << "do_fill: " << +do_fill;
										cout << " lower level wq is full!" << " fill_addr: " << hex << WQ.entry[index].address;
										cout << " victim_addr: " << block[set][way].tag << dec << endl; });
							}
							else { 
								PACKET writeback_packet;

								writeback_packet.fill_level = fill_level << 1;
								writeback_packet.cpu = writeback_cpu;
								writeback_packet.address = block[set][way].address;
								writeback_packet.full_addr = block[set][way].full_addr;
								writeback_packet.data = block[set][way].data;
								writeback_packet.instr_id = WQ.entry[index].instr_id;
								writeback_packet.ip = 0;
								writeback_packet.type = WRITEBACK;
								writeback_packet.event_cycle = current_core_cycle[writeback_cpu];

								lower_level->add_wq(&writeback_packet);
							}
						}
#ifdef SANITY_CHECK
						else {
							// sanity check
							if (cache_type != IS_STLB)
								assert(0);
						}
#endif
					}

					if (do_fill) {
						// update prefetcher
						if (cache_type == IS_L1I)
							l1i_prefetcher_cache_fill(writeback_cpu, ((WQ.entry[index].ip)>>LOG2_BLOCK_SIZE)<<LOG2_BLOCK_SIZE, set, way, 0, ((block[set][way].ip)>>LOG2_BLOCK_SIZE)<<LOG2_BLOCK_SIZE);
						if (cache_type == IS_L1D)
							l1d_prefetcher_cache_fill(WQ.entry[index].full_addr, set, way, 0, block[set][way].address<<LOG2_BLOCK_SIZE, WQ.entry[index].pf_metadata);
						else if (cache_type == IS_L2C)
							WQ.entry[index].pf_metadata = l2c_prefetcher_cache_fill(WQ.entry[index].address<<LOG2_BLOCK_SIZE, set, way, 0,
									block[set][way].address<<LOG2_BLOCK_SIZE, WQ.entry[index].pf_metadata);
						if (cache_type == IS_LLC)
						{
							cpu = writeback_cpu;
							WQ.entry[index].pf_metadata =llc_prefetcher_cache_fill(WQ.entry[index].address<<LOG2_BLOCK_SIZE, set, way, 0,
									block[set][way].address<<LOG2_BLOCK_SIZE, WQ.entry[index].pf_metadata);
							cpu = 0;
						}

						// update replacement policy
						if (cache_type == IS_LLC) {
							llc_update_replacement_state(writeback_cpu, set, way, WQ.entry[index].full_addr, WQ.entry[index].ip, block[set][way].full_addr, WQ.entry[index].type, 0);
						}
						else
							update_replacement_state(writeback_cpu, set, way, WQ.entry[index].full_addr, WQ.entry[index].ip, block[set][way].full_addr, WQ.entry[index].type, 0);

						// COLLECT STATS
						sim_miss[writeback_cpu][WQ.entry[index].type]++;
						sim_access[writeback_cpu][WQ.entry[index].type]++;

						fill_cache(set, way, &WQ.entry[index]);

						// mark dirty
						block[set][way].dirty = 1; 

						// check fill level
						if (WQ.entry[index].fill_level < fill_level) {

							if(fill_level == FILL_L2)
							{
								if(WQ.entry[index].fill_l1i)
								{
									upper_level_icache[writeback_cpu]->return_data(&WQ.entry[index]);
								}
								if(WQ.entry[index].fill_l1d)
								{
									upper_level_dcache[writeback_cpu]->return_data(&WQ.entry[index]);
								}
							}
							else
							{
								if (WQ.entry[index].instruction)
									upper_level_icache[writeback_cpu]->return_data(&WQ.entry[index]);
								if (WQ.entry[index].is_data)
									upper_level_dcache[writeback_cpu]->return_data(&WQ.entry[index]);
							}
						}

						MISS[WQ.entry[index].type]++;
						ACCESS[WQ.entry[index].type]++;

						// remove this entry from WQ
						WQ.remove_queue(&WQ.entry[index]);
					}
				}
			}
		}
	}

	void CACHE::handle_read()
	{
		// handle read
		for (uint32_t i=0; i<MAX_READ; i++) {

			uint32_t read_cpu = RQ.entry[RQ.head].cpu;
			if (read_cpu == NUM_CPUS)
				return;

			// handle the oldest entry
			if ((RQ.entry[RQ.head].event_cycle <= current_core_cycle[read_cpu]) && (RQ.occupancy > 0)) {
				int index = RQ.head;

				// access cache
				uint32_t set = get_set(RQ.entry[index].address);
				int way = check_hit(&RQ.entry[index]);

				//if (cache_type == IS_ITLB)
				//	ooo_cpu[0].STLB.code_footprint.insert(make_pair(RQ.entry[index].address, 1));


				if (way >= 0) { // read hit

					if (cache_type == IS_ITLB) {
						RQ.entry[index].instruction_pa = block[set][way].data;
						if (PROCESSED.occupancy < PROCESSED.SIZE)
							PROCESSED.add_queue(&RQ.entry[index]);
					}
					else if (cache_type == IS_DTLB) {
						RQ.entry[index].data_pa = block[set][way].data;
						if (PROCESSED.occupancy < PROCESSED.SIZE)
							PROCESSED.add_queue(&RQ.entry[index]);
					}
					else if (cache_type == IS_STLB){
						RQ.entry[index].data = block[set][way].data;

						if(current_core_cycle[read_cpu] < block[set][way].stalls)
							add_stall_prefetch(block[set][way].stalls-current_core_cycle[read_cpu], read_cpu);
						int bits, rowhit=-1, victim, iflag = 0;
						

						pair<int, int> answer = make_pair(-1,-1);

						int * free_indexes;
						free_indexes = sorted_free_distances();

						if(unsigned(RQ.entry[index].instruction) != 0)
							iflag=1;

						if(iflag == 0){ //	Was iflag==1
							free_indexes = sorted_free_distances();
							stlb_prefetcher_operate(RQ.entry[index].address, RQ.entry[index].ip, 0, RQ.entry[index].type, answer.first, warmup_complete[cpu], free_indexes, RQ.entry[index].instr_id, 1, current_core_cycle[read_cpu]);
							stlb_prefetcher_cache_fill(RQ.entry[index].address, 0, 0, 0, 0);
						}
					}
					else if (cache_type == IS_L1I) {
						if (PROCESSED.occupancy < PROCESSED.SIZE)
							PROCESSED.add_queue(&RQ.entry[index]);
					}
					//else if (cache_type == IS_L1D)
					else if ((cache_type == IS_L1D) && (RQ.entry[index].type != PREFETCH)) {
						if (PROCESSED.occupancy < PROCESSED.SIZE)
							PROCESSED.add_queue(&RQ.entry[index]);
					}

					// update prefetcher on load instruction
					if (RQ.entry[index].type == LOAD) {
						if(cache_type == IS_L1I)
							l1i_prefetcher_cache_operate(read_cpu, RQ.entry[index].ip, 1, block[set][way].prefetch);
						if (cache_type == IS_L1D) 
							l1d_prefetcher_operate(RQ.entry[index].full_addr, RQ.entry[index].ip, 1, RQ.entry[index].type);
						else if (cache_type == IS_L2C)
							l2c_prefetcher_operate(block[set][way].address<<LOG2_BLOCK_SIZE, RQ.entry[index].ip, 1, RQ.entry[index].type, 0);
						else if (cache_type == IS_LLC)
						{
							cpu = read_cpu;
							llc_prefetcher_operate(block[set][way].address<<LOG2_BLOCK_SIZE, RQ.entry[index].ip, 1, RQ.entry[index].type, 0);
							cpu = 0;
						}
					}

					// update replacement policy
					if (cache_type == IS_LLC) {
						llc_update_replacement_state(read_cpu, set, way, block[set][way].full_addr, RQ.entry[index].ip, 0, RQ.entry[index].type, 1);

					}
					else
						update_replacement_state(read_cpu, set, way, block[set][way].full_addr, RQ.entry[index].ip, 0, RQ.entry[index].type, 1);

					// COLLECT STATS
					sim_hit[read_cpu][RQ.entry[index].type]++;
					sim_access[read_cpu][RQ.entry[index].type]++;

					// check fill level
					if (RQ.entry[index].fill_level < fill_level) {

						if(fill_level == FILL_L2)
						{
							if(RQ.entry[index].fill_l1i)
							{
								upper_level_icache[read_cpu]->return_data(&RQ.entry[index]);
							}
							if(RQ.entry[index].fill_l1d)
							{
								upper_level_dcache[read_cpu]->return_data(&RQ.entry[index]);
							}
						}
						else
						{
							if (RQ.entry[index].instruction)
								upper_level_icache[read_cpu]->return_data(&RQ.entry[index]);
							if (RQ.entry[index].is_data)
								upper_level_dcache[read_cpu]->return_data(&RQ.entry[index]);
						}
					}

					// update prefetch stats and reset prefetch bit
					if (block[set][way].prefetch) {
						pf_useful++;
						block[set][way].prefetch = 0;
					}
					block[set][way].used = 1;

					HIT[RQ.entry[index].type]++;
					ACCESS[RQ.entry[index].type]++;

					// remove this entry from RQ
					RQ.remove_queue(&RQ.entry[index]);
					reads_available_this_cycle--;
				}
				else { // read miss

					DP ( if (warmup_complete[read_cpu]) {
							cout << "[" << NAME << "] " << __func__ << " read miss";
							cout << " instr_id: " << RQ.entry[index].instr_id << " address: " << hex << RQ.entry[index].address;
							cout << " full_addr: " << RQ.entry[index].full_addr << dec;
							cout << " cycle: " << RQ.entry[index].event_cycle << endl; });

					// check mshr
					uint8_t miss_handled = 1;
					int mshr_index = check_mshr(&RQ.entry[index]);

					if(mshr_index == -2)
					{
						// this is a data/instruction collision in the MSHR, so we have to wait before we can allocate this miss
						miss_handled = 0;
					}
					else if ((mshr_index == -1) && (MSHR.occupancy < MSHR_SIZE)) { // this is a new miss

						if(cache_type == IS_LLC)
						{
							// check to make sure the DRAM RQ has room for this LLC read miss
							if (lower_level->get_occupancy(1, RQ.entry[index].address) == lower_level->get_size(1, RQ.entry[index].address))
							{
								miss_handled = 0;
							}
							else
							{
								add_mshr(&RQ.entry[index]);
								if(lower_level)
								{
									lower_level->add_rq(&RQ.entry[index]);
								}
							}
						}
						else
						{
							add_mshr(&RQ.entry[index]);

							if (lower_level)
								lower_level->add_rq(&RQ.entry[index]);
							else {
								if (cache_type == IS_STLB) {
									uint64_t pa, current_vpn = RQ.entry[index].address;
									int bits, rowhit=-1, victim, iflag = 0;
									pair<int, int> answer;

									int * free_indexes;
									free_indexes = sorted_free_distances();

									if(unsigned(RQ.entry[index].instruction) != 0){
										stlb_misses[1]++;
										iflag = 1;
										decay_timer++;
										decay_conf_timer++;
									}
									else{
										stlb_misses[0]++;
									}

									int fctb_found_pos = -10;
									if(iflag == 1){
										refresh_fctb(current_core_cycle[read_cpu]);
										fctb_found_pos = search_fctb(RQ.entry[index].address);
										answer = check_hit_stlb_pq(RQ.entry[index].address);
									}
									else{
										answer = make_pair(-1,-1);
									}

									pair<uint64_t, uint64_t> v2p;
									if(answer.first == -1){
										if(iflag){
											if(fctb_found_pos == -10){
												v2p = va_to_pa(read_cpu, RQ.entry[index].instr_id, RQ.entry[index].full_addr, RQ.entry[index].address, RQ.entry[index].ip, RQ.entry[index].type, iflag, 0);
												pa = v2p.first;
												int victim_entry = fctb_replacement_policy();
												fctb[victim_entry][0] = current_vpn;
												fctb[victim_entry][1] = (RQ.entry[index].full_addr & 0x7000)/4096;
												fctb[victim_entry][2] = current_core_cycle[read_cpu];
												fctb[victim_entry][3] = v2p.second;

											}
											else{
												pa = va_to_pa_prefetch(read_cpu, RQ.entry[index].full_addr, RQ.entry[index].address);
												if(pa == 0){
													v2p = va_to_pa(read_cpu, RQ.entry[index].instr_id, RQ.entry[index].full_addr, RQ.entry[index].address, RQ.entry[index].ip, RQ.entry[index].type, iflag, 0);
													pa = v2p.first;
												}
											}
										}
										else{
											v2p = va_to_pa(read_cpu, RQ.entry[index].instr_id, RQ.entry[index].full_addr, RQ.entry[index].address, RQ.entry[index].ip, RQ.entry[index].type, iflag, 0);
											pa  = v2p.first;
										}
										if (iflag == 1)
											pf_misses_pq++;
									}
									else{
										if(PQ.entry[answer.first].free_bit == 1 && PQ.entry[answer.first].free_distance != 0){
											rfhits[1]++;
											free_hits[PQ.entry[answer.first].free_distance + 6 + (PQ.entry[answer.first].free_distance < 0)*1]++;
										}
										else
											rfhits[0]++;

										if (warmup_complete[cpu]){
											uint64_t num_cycles;

											if(PQ.entry[answer.first].irip)
												irip_hits++;
											else
												sdp_hits++;

											if(answer.second == 0){
												num_cycles = current_core_cycle[read_cpu] - PQ.return_event_cycle(answer.first);
												if (num_cycles < PQ.entry[answer.first].stall_cycles){
													add_stall_prefetch(PQ.entry[answer.first].stall_cycles - num_cycles, read_cpu);
												}
											}
											else{
												cout << "you should not be here since we are using only one PQ" << endl;
											}
										}

										if(answer.second == 0){
											pa = PQ.entry[answer.first].data;
											if(PQ.entry[answer.first].lad == 1)
												hit_prefetches_lad++;
											PQ.remove_queue_prefetch(&PQ.entry[answer.first], answer.first);
										}
										else{
											cout << "you should not be here since we are using only one PQ" << endl;
										}

										pf_hits_pq++;

									}

									RQ.entry[index].data = pa >> LOG2_PAGE_SIZE; 
									RQ.entry[index].event_cycle = current_core_cycle[read_cpu];
									return_data(&RQ.entry[index]);

									if(iflag == 0){ //changed to dSTLB
										free_indexes = sorted_free_distances();
										stlb_prefetcher_operate(RQ.entry[index].address, RQ.entry[index].ip, 0, RQ.entry[index].type, answer.first, warmup_complete[cpu], free_indexes, RQ.entry[index].instr_id, 0, current_core_cycle[read_cpu]);
										stlb_prefetcher_cache_fill(RQ.entry[index].address, 0, 0, 0, 0);
									}
								}
							}
						}
					}
					else {
						if ((mshr_index == -1) && (MSHR.occupancy == MSHR_SIZE)) { // not enough MSHR resource

							// cannot handle miss request until one of MSHRs is available
							miss_handled = 0;
							STALL[RQ.entry[index].type]++;
						}
						else if (mshr_index != -1) { // already in-flight miss

							// mark merged consumer
							if (RQ.entry[index].type == RFO) {

								if (RQ.entry[index].tlb_access) {
									uint32_t sq_index = RQ.entry[index].sq_index;
									MSHR.entry[mshr_index].store_merged = 1;
									MSHR.entry[mshr_index].sq_index_depend_on_me.insert (sq_index);
									MSHR.entry[mshr_index].sq_index_depend_on_me.join (RQ.entry[index].sq_index_depend_on_me, SQ_SIZE);
								}

								if (RQ.entry[index].load_merged) {
									//uint32_t lq_index = RQ.entry[index].lq_index; 
									MSHR.entry[mshr_index].load_merged = 1;
									//MSHR.entry[mshr_index].lq_index_depend_on_me[lq_index] = 1;
									MSHR.entry[mshr_index].lq_index_depend_on_me.join (RQ.entry[index].lq_index_depend_on_me, LQ_SIZE);
								}
							}
							else {
								if (RQ.entry[index].instruction) {
									uint32_t rob_index = RQ.entry[index].rob_index;
									MSHR.entry[mshr_index].instruction = 1; // add as instruction type
									MSHR.entry[mshr_index].instr_merged = 1;
									MSHR.entry[mshr_index].rob_index_depend_on_me.insert (rob_index);

									DP (if (warmup_complete[MSHR.entry[mshr_index].cpu]) {
											cout << "[INSTR_MERGED] " << __func__ << " cpu: " << MSHR.entry[mshr_index].cpu << " instr_id: " << MSHR.entry[mshr_index].instr_id;
											cout << " merged rob_index: " << rob_index << " instr_id: " << RQ.entry[index].instr_id << endl; });

									if (RQ.entry[index].instr_merged) {
										MSHR.entry[mshr_index].rob_index_depend_on_me.join (RQ.entry[index].rob_index_depend_on_me, ROB_SIZE);
										DP (if (warmup_complete[MSHR.entry[mshr_index].cpu]) {
												cout << "[INSTR_MERGED] " << __func__ << " cpu: " << MSHR.entry[mshr_index].cpu << " instr_id: " << MSHR.entry[mshr_index].instr_id;
												cout << " merged rob_index: " << i << " instr_id: N/A" << endl; });
									}
								}
								else 
								{
									uint32_t lq_index = RQ.entry[index].lq_index;
									MSHR.entry[mshr_index].is_data = 1; // add as data type
									MSHR.entry[mshr_index].load_merged = 1;
									MSHR.entry[mshr_index].lq_index_depend_on_me.insert (lq_index);

									DP (if (warmup_complete[read_cpu]) {
											cout << "[DATA_MERGED] " << __func__ << " cpu: " << read_cpu << " instr_id: " << RQ.entry[index].instr_id;
											cout << " merged rob_index: " << RQ.entry[index].rob_index << " instr_id: " << RQ.entry[index].instr_id << " lq_index: " << RQ.entry[index].lq_index << endl; });
									MSHR.entry[mshr_index].lq_index_depend_on_me.join (RQ.entry[index].lq_index_depend_on_me, LQ_SIZE);
									if (RQ.entry[index].store_merged) {
										MSHR.entry[mshr_index].store_merged = 1;
										MSHR.entry[mshr_index].sq_index_depend_on_me.join (RQ.entry[index].sq_index_depend_on_me, SQ_SIZE);
									}
								}
							}

							// update fill_level
							if (RQ.entry[index].fill_level < MSHR.entry[mshr_index].fill_level)
								MSHR.entry[mshr_index].fill_level = RQ.entry[index].fill_level;

							if((RQ.entry[index].fill_l1i) && (MSHR.entry[mshr_index].fill_l1i != 1))
							{
								MSHR.entry[mshr_index].fill_l1i = 1;
							}
							if((RQ.entry[index].fill_l1d) && (MSHR.entry[mshr_index].fill_l1d != 1))
							{
								MSHR.entry[mshr_index].fill_l1d = 1;
							}

							// update request
							if (MSHR.entry[mshr_index].type == PREFETCH) {
								uint8_t  prior_returned = MSHR.entry[mshr_index].returned;
								uint64_t prior_event_cycle = MSHR.entry[mshr_index].event_cycle;
								MSHR.entry[mshr_index] = RQ.entry[index];

								// in case request is already returned, we should keep event_cycle and retunred variables
								MSHR.entry[mshr_index].returned = prior_returned;
								MSHR.entry[mshr_index].event_cycle = prior_event_cycle;
							}

							MSHR_MERGED[RQ.entry[index].type]++;

							DP ( if (warmup_complete[read_cpu]) {
									cout << "[" << NAME << "] " << __func__ << " mshr merged";
									cout << " instr_id: " << RQ.entry[index].instr_id << " prior_id: " << MSHR.entry[mshr_index].instr_id; 
									cout << " address: " << hex << RQ.entry[index].address;
									cout << " full_addr: " << RQ.entry[index].full_addr << dec;
									cout << " cycle: " << RQ.entry[index].event_cycle << endl; });
						}
						else { // WE SHOULD NOT REACH HERE
							cerr << "[" << NAME << "] MSHR errors" << endl;
							assert(0);
						}
					}

					if (miss_handled) {
						// update prefetcher on load instruction
						if (RQ.entry[index].type == LOAD) {
							if(cache_type == IS_L1I)
								l1i_prefetcher_cache_operate(read_cpu, RQ.entry[index].ip, 0, 0);
							if (cache_type == IS_L1D) 
								l1d_prefetcher_operate(RQ.entry[index].full_addr, RQ.entry[index].ip, 0, RQ.entry[index].type);
							if (cache_type == IS_L2C)
								l2c_prefetcher_operate(RQ.entry[index].address<<LOG2_BLOCK_SIZE, RQ.entry[index].ip, 0, RQ.entry[index].type, 0);
							if (cache_type == IS_LLC)
							{
								cpu = read_cpu;
								llc_prefetcher_operate(RQ.entry[index].address<<LOG2_BLOCK_SIZE, RQ.entry[index].ip, 0, RQ.entry[index].type, 0);
								cpu = 0;
							}
						}

						MISS[RQ.entry[index].type]++;
						ACCESS[RQ.entry[index].type]++;

						// remove this entry from RQ
						RQ.remove_queue(&RQ.entry[index]);
						reads_available_this_cycle--;
					}
				}
			}
			else
			{
				return;
			}

			if(reads_available_this_cycle == 0)
			{
				return;
			}
		}
	}

	void CACHE::handle_prefetch()
	{
		// handle prefetch

		for (uint32_t i=0; i<MAX_READ; i++) {

			uint32_t prefetch_cpu = PQ.entry[PQ.head].cpu;
			if (prefetch_cpu == NUM_CPUS)
				return;

			// handle the oldest entry
			if ((PQ.entry[PQ.head].event_cycle <= current_core_cycle[prefetch_cpu]) && (PQ.occupancy > 0)) {
				int index = PQ.head;

				// access cache
				uint32_t set = get_set(PQ.entry[index].address);
				int way = check_hit(&PQ.entry[index]);

				if (way >= 0) { // prefetch hit

					// update replacement policy
					if (cache_type == IS_LLC) {
						llc_update_replacement_state(prefetch_cpu, set, way, block[set][way].full_addr, PQ.entry[index].ip, 0, PQ.entry[index].type, 1);

					}
					else
						update_replacement_state(prefetch_cpu, set, way, block[set][way].full_addr, PQ.entry[index].ip, 0, PQ.entry[index].type, 1);

					// COLLECT STATS
					sim_hit[prefetch_cpu][PQ.entry[index].type]++;
					sim_access[prefetch_cpu][PQ.entry[index].type]++;

					// run prefetcher on prefetches from higher caches
					if(PQ.entry[index].pf_origin_level < fill_level)
					{
						if (cache_type == IS_L1D)
							l1d_prefetcher_operate(PQ.entry[index].full_addr, PQ.entry[index].ip, 1, PREFETCH);
						else if (cache_type == IS_L2C)
							PQ.entry[index].pf_metadata = l2c_prefetcher_operate(block[set][way].address<<LOG2_BLOCK_SIZE, PQ.entry[index].ip, 1, PREFETCH, PQ.entry[index].pf_metadata);
						else if (cache_type == IS_LLC)
						{
							cpu = prefetch_cpu;
							PQ.entry[index].pf_metadata = llc_prefetcher_operate(block[set][way].address<<LOG2_BLOCK_SIZE, PQ.entry[index].ip, 1, PREFETCH, PQ.entry[index].pf_metadata);
							cpu = 0;
						}
					}

					// check fill level
					if (PQ.entry[index].fill_level < fill_level) {

						if(fill_level == FILL_L2)
						{
							if(PQ.entry[index].fill_l1i)
							{
								upper_level_icache[prefetch_cpu]->return_data(&PQ.entry[index]);
							}
							if(PQ.entry[index].fill_l1d)
							{
								upper_level_dcache[prefetch_cpu]->return_data(&PQ.entry[index]);
							}
						}
						else
						{
							if (PQ.entry[index].instruction)
								upper_level_icache[prefetch_cpu]->return_data(&PQ.entry[index]);
							if (PQ.entry[index].is_data)
								upper_level_dcache[prefetch_cpu]->return_data(&PQ.entry[index]);
						}
					}

					HIT[PQ.entry[index].type]++;
					ACCESS[PQ.entry[index].type]++;

					// remove this entry from PQ
					PQ.remove_queue(&PQ.entry[index]);
					reads_available_this_cycle--;
				}
				else { // prefetch miss

					DP ( if (warmup_complete[prefetch_cpu]) {
							cout << "[" << NAME << "] " << __func__ << " prefetch miss";
							cout << " instr_id: " << PQ.entry[index].instr_id << " address: " << hex << PQ.entry[index].address;
							cout << " full_addr: " << PQ.entry[index].full_addr << dec << " fill_level: " << PQ.entry[index].fill_level;
							cout << " cycle: " << PQ.entry[index].event_cycle << endl; });

					// check mshr
					uint8_t miss_handled = 1;
					int mshr_index = check_mshr(&PQ.entry[index]);

					if(mshr_index == -2)
					{
						// this is a data/instruction collision in the MSHR, so we have to wait before we can allocate this miss
						miss_handled = 0;
					}
					else if ((mshr_index == -1) && (MSHR.occupancy < MSHR_SIZE)) { // this is a new miss

						DP ( if (warmup_complete[PQ.entry[index].cpu]) {
								cout << "[" << NAME << "_PQ] " <<  __func__ << " want to add instr_id: " << PQ.entry[index].instr_id << " address: " << hex << PQ.entry[index].address;
								cout << " full_addr: " << PQ.entry[index].full_addr << dec;
								cout << " occupancy: " << lower_level->get_occupancy(3, PQ.entry[index].address) << " SIZE: " << lower_level->get_size(3, PQ.entry[index].address) << endl; });

						// first check if the lower level PQ is full or not
						// this is possible since multiple prefetchers can exist at each level of caches
						if (lower_level) {
							if (cache_type == IS_LLC) {
								if (lower_level->get_occupancy(1, PQ.entry[index].address) == lower_level->get_size(1, PQ.entry[index].address))
									miss_handled = 0;
								else {

									// run prefetcher on prefetches from higher caches
									if(PQ.entry[index].pf_origin_level < fill_level)
									{
										if (cache_type == IS_LLC)
										{
											cpu = prefetch_cpu;
											PQ.entry[index].pf_metadata = llc_prefetcher_operate(PQ.entry[index].address<<LOG2_BLOCK_SIZE, PQ.entry[index].ip, 0, PREFETCH, PQ.entry[index].pf_metadata);
											cpu = 0;
										}
									}

									// add it to MSHRs if this prefetch miss will be filled to this cache level
									if (PQ.entry[index].fill_level <= fill_level)
										add_mshr(&PQ.entry[index]);

									lower_level->add_rq(&PQ.entry[index]); // add it to the DRAM RQ
								}
							}
							else {
								if (lower_level->get_occupancy(3, PQ.entry[index].address) == lower_level->get_size(3, PQ.entry[index].address))
									miss_handled = 0;
								else {

									// run prefetcher on prefetches from higher caches
									if(PQ.entry[index].pf_origin_level < fill_level)
									{
										if (cache_type == IS_L1D)
											l1d_prefetcher_operate(PQ.entry[index].full_addr, PQ.entry[index].ip, 0, PREFETCH);
										if (cache_type == IS_L2C)
											PQ.entry[index].pf_metadata = l2c_prefetcher_operate(PQ.entry[index].address<<LOG2_BLOCK_SIZE, PQ.entry[index].ip, 0, PREFETCH, PQ.entry[index].pf_metadata);
									}

									// add it to MSHRs if this prefetch miss will be filled to this cache level
									if (PQ.entry[index].fill_level <= fill_level)
										add_mshr(&PQ.entry[index]);

									lower_level->add_pq(&PQ.entry[index]); // add it to the DRAM RQ
								}
							}
						}
					}
					else {
						if ((mshr_index == -1) && (MSHR.occupancy == MSHR_SIZE)) { // not enough MSHR resource

							// TODO: should we allow prefetching with lower fill level at this case?

							// cannot handle miss request until one of MSHRs is available
							miss_handled = 0;
							STALL[PQ.entry[index].type]++;
						}
						else if (mshr_index != -1) { // already in-flight miss

							// no need to update request except fill_level
							// update fill_level
							if (PQ.entry[index].fill_level < MSHR.entry[mshr_index].fill_level)
								MSHR.entry[mshr_index].fill_level = PQ.entry[index].fill_level;

							if((PQ.entry[index].fill_l1i) && (MSHR.entry[mshr_index].fill_l1i != 1))
							{
								MSHR.entry[mshr_index].fill_l1i = 1;
							}
							if((PQ.entry[index].fill_l1d) && (MSHR.entry[mshr_index].fill_l1d != 1))
							{
								MSHR.entry[mshr_index].fill_l1d = 1;
							}

							MSHR_MERGED[PQ.entry[index].type]++;

							DP ( if (warmup_complete[prefetch_cpu]) {
									cout << "[" << NAME << "] " << __func__ << " mshr merged";
									cout << " instr_id: " << PQ.entry[index].instr_id << " prior_id: " << MSHR.entry[mshr_index].instr_id; 
									cout << " address: " << hex << PQ.entry[index].address;
									cout << " full_addr: " << PQ.entry[index].full_addr << dec << " fill_level: " << MSHR.entry[mshr_index].fill_level;
									cout << " cycle: " << MSHR.entry[mshr_index].event_cycle << endl; });
						}
						else { // WE SHOULD NOT REACH HERE
							cerr << "[" << NAME << "] MSHR errors" << endl;
							assert(0);
						}
					}

					if (miss_handled) {

						DP ( if (warmup_complete[prefetch_cpu]) {
								cout << "[" << NAME << "] " << __func__ << " prefetch miss handled";
								cout << " instr_id: " << PQ.entry[index].instr_id << " address: " << hex << PQ.entry[index].address;
								cout << " full_addr: " << PQ.entry[index].full_addr << dec << " fill_level: " << PQ.entry[index].fill_level;
								cout << " cycle: " << PQ.entry[index].event_cycle << endl; });

						MISS[PQ.entry[index].type]++;
						ACCESS[PQ.entry[index].type]++;

						// remove this entry from PQ
						PQ.remove_queue(&PQ.entry[index]);
						reads_available_this_cycle--;
					}
				}
			}
			else
			{
				return;
			}

			if(reads_available_this_cycle == 0)
			{
				return;
			}
		}
	}

	void CACHE::operate()
	{
		handle_fill();
		handle_writeback();
		reads_available_this_cycle = MAX_READ;
		handle_read();

		if (PQ.occupancy && (reads_available_this_cycle > 0))
			if(cache_type != IS_STLB)
				handle_prefetch();
	}

	uint32_t CACHE::get_set(uint64_t address)
	{
		return (uint32_t) (address & ((1 << lg2(NUM_SET)) - 1)); 
	}

	uint32_t CACHE::get_way(uint64_t address, uint32_t set)
	{
		for (uint32_t way=0; way<NUM_WAY; way++) {
			if (block[set][way].valid && (block[set][way].tag == address)) 
				return way;
		}

		return NUM_WAY;
	}

	void CACHE::fill_cache(uint32_t set, uint32_t way, PACKET *packet)
	{
#ifdef SANITY_CHECK
		if (cache_type == IS_ITLB) {
			if (packet->data == 0){
				assert(0);
			}
		}

		if (cache_type == IS_DTLB) {
			if (packet->data == 0){
				assert(0);
			}
		}

		if (cache_type == IS_STLB) {
			if (packet->data == 0){
				assert(0);
			}
		}
#endif
		if (block[set][way].prefetch && (block[set][way].used == 0))
			pf_useless++;

		if (block[set][way].valid == 0)
			block[set][way].valid = 1;
		block[set][way].dirty = 0;
		block[set][way].prefetch = (packet->type == PREFETCH) ? 1 : 0;
		block[set][way].used = 0;

		if (block[set][way].prefetch)
			pf_fill++;

		block[set][way].delta = packet->delta;
		block[set][way].depth = packet->depth;
		block[set][way].signature = packet->signature;
		block[set][way].confidence = packet->confidence;

		block[set][way].tag = packet->address;
		block[set][way].address = packet->address;
		block[set][way].full_addr = packet->full_addr;
		block[set][way].data = packet->data;
		block[set][way].ip = packet->ip;
		block[set][way].cpu = packet->cpu;
		block[set][way].instr_id = packet->instr_id;

		DP ( if (warmup_complete[packet->cpu]) {
				cout << "[" << NAME << "] " << __func__ << " set: " << set << " way: " << way;
				cout << " lru: " << block[set][way].lru << " tag: " << hex << block[set][way].tag << " full_addr: " << block[set][way].full_addr;
				cout << " data: " << block[set][way].data << dec << endl; });
	}

	pair<int, int> CACHE::check_hit_stlb_pq(uint64_t vpn){
		/* In case of hit in a PQ, we return a tuple with the position of the found translation and the identifier of the corresponding PQ */
		pair<int, uint64_t> temp_pair  = PQ.check_queue_vpn(vpn);
		pair<int, int> return_pair;

		if(temp_pair.first == -1){
			return_pair = make_pair(-1, -1);
		}
		else
			return_pair = make_pair(temp_pair.first,0);

		return return_pair;
	}



	int CACHE::check_hit(PACKET *packet)
	{
		uint32_t set = get_set(packet->address);
		int match_way = -1;

		if (NUM_SET < set) {
			cerr << "[" << NAME << "_ERROR] " << __func__ << " invalid set index: " << set << " NUM_SET: " << NUM_SET;
			cerr << " address: " << hex << packet->address << " full_addr: " << packet->full_addr << dec;
			cerr << " event: " << packet->event_cycle << endl;
			assert(0);
		}

		// hit
		for (uint32_t way=0; way<NUM_WAY; way++) {
			if (block[set][way].valid && (block[set][way].tag == packet->address)) {

				match_way = way;

				DP ( if (warmup_complete[packet->cpu]) {
						cout << "[" << NAME << "] " << __func__ << " instr_id: " << packet->instr_id << " type: " << +packet->type << hex << " addr: " << packet->address;
						cout << " full_addr: " << packet->full_addr << " tag: " << block[set][way].tag << " data: " << block[set][way].data << dec;
						cout << " set: " << set << " way: " << way << " lru: " << block[set][way].lru;
						cout << " event: " << packet->event_cycle << " cycle: " << current_core_cycle[cpu] << endl; });

				break;
			}
		}

		return match_way;
	}

	int CACHE::invalidate_entry(uint64_t inval_addr)
	{
		uint32_t set = get_set(inval_addr);
		int match_way = -1;

		if (NUM_SET < set) {
			cerr << "[" << NAME << "_ERROR] " << __func__ << " invalid set index: " << set << " NUM_SET: " << NUM_SET;
			cerr << " inval_addr: " << hex << inval_addr << dec << endl;
			assert(0);
		}

		// invalidate
		for (uint32_t way=0; way<NUM_WAY; way++) {
			if (block[set][way].valid && (block[set][way].tag == inval_addr)) {

				block[set][way].valid = 0;

				match_way = way;

				DP ( if (warmup_complete[cpu]) {
						cout << "[" << NAME << "] " << __func__ << " inval_addr: " << hex << inval_addr;  
						cout << " tag: " << block[set][way].tag << " data: " << block[set][way].data << dec;
						cout << " set: " << set << " way: " << way << " lru: " << block[set][way].lru << " cycle: " << current_core_cycle[cpu] << endl; });

				break;
			}
		}

		return match_way;
	}

	int CACHE::add_rq(PACKET *packet)
	{
		// check for the latest wirtebacks in the write queue
		int wq_index = WQ.check_queue(packet);
		if (wq_index != -1) {

			// check fill level
			if (packet->fill_level < fill_level) {

				packet->data = WQ.entry[wq_index].data;

				if(fill_level == FILL_L2)
				{
					if(packet->fill_l1i)
					{
						upper_level_icache[packet->cpu]->return_data(packet);
					}
					if(packet->fill_l1d)
					{
						upper_level_dcache[packet->cpu]->return_data(packet);
					}
				}
				else
				{
					if (packet->instruction)
						upper_level_icache[packet->cpu]->return_data(packet);
					if (packet->is_data)
						upper_level_dcache[packet->cpu]->return_data(packet);
				}
			}

#ifdef SANITY_CHECK
			if (cache_type == IS_ITLB)
				assert(0);
			else if (cache_type == IS_DTLB)
				assert(0);
			else if (cache_type == IS_L1I)
				assert(0);
#endif
			// update processed packets
			if ((cache_type == IS_L1D) && (packet->type != PREFETCH)) {
				if (PROCESSED.occupancy < PROCESSED.SIZE)
					PROCESSED.add_queue(packet);

				DP ( if (warmup_complete[packet->cpu]) {
						cout << "[" << NAME << "_RQ] " << __func__ << " instr_id: " << packet->instr_id << " found recent writebacks";
						cout << hex << " read: " << packet->address << " writeback: " << WQ.entry[wq_index].address << dec;
						cout << " index: " << MAX_READ << " rob_signal: " << packet->rob_signal << endl; });
			}

			HIT[packet->type]++;
			ACCESS[packet->type]++;

			WQ.FORWARD++;
			RQ.ACCESS++;

			return -1;
		}

		// check for duplicates in the read queue
		int index = RQ.check_queue(packet);
		if (index != -1) {

			if (packet->instruction) {
				uint32_t rob_index = packet->rob_index;
				RQ.entry[index].rob_index_depend_on_me.insert (rob_index);
				RQ.entry[index].instruction = 1; // add as instruction type
				RQ.entry[index].instr_merged = 1;

				DP (if (warmup_complete[packet->cpu]) {
						cout << "[INSTR_MERGED] " << __func__ << " cpu: " << packet->cpu << " instr_id: " << RQ.entry[index].instr_id;
						cout << " merged rob_index: " << rob_index << " instr_id: " << packet->instr_id << endl; });
			}
			else 
			{
				// mark merged consumer
				if (packet->type == RFO) {

					uint32_t sq_index = packet->sq_index;
					RQ.entry[index].sq_index_depend_on_me.insert (sq_index);
					RQ.entry[index].store_merged = 1;
				}
				else {
					uint32_t lq_index = packet->lq_index; 
					RQ.entry[index].lq_index_depend_on_me.insert (lq_index);
					RQ.entry[index].load_merged = 1;

					DP (if (warmup_complete[packet->cpu]) {
							cout << "[DATA_MERGED] " << __func__ << " cpu: " << packet->cpu << " instr_id: " << RQ.entry[index].instr_id;
							cout << " merged rob_index: " << packet->rob_index << " instr_id: " << packet->instr_id << " lq_index: " << packet->lq_index << endl; });
				}
				RQ.entry[index].is_data = 1; // add as data type
			}

			if((packet->fill_l1i) && (RQ.entry[index].fill_l1i != 1))
			{
				RQ.entry[index].fill_l1i = 1;
			}
			if((packet->fill_l1d) && (RQ.entry[index].fill_l1d != 1))
			{
				RQ.entry[index].fill_l1d = 1;
			}

			RQ.MERGED++;
			RQ.ACCESS++;

			return index; // merged index
		}

		// check occupancy
		if (RQ.occupancy == RQ_SIZE) {
			RQ.FULL++;

			return -2; // cannot handle this request
		}

		// if there is no duplicate, add it to RQ
		index = RQ.tail;

#ifdef SANITY_CHECK
		if (RQ.entry[index].address != 0) {
			cerr << "[" << NAME << "_ERROR] " << __func__ << " is not empty index: " << index;
			cerr << " address: " << hex << RQ.entry[index].address;
			cerr << " full_addr: " << RQ.entry[index].full_addr << dec << endl;
			assert(0);
		}
#endif

		RQ.entry[index] = *packet;

		// ADD LATENCY
		if (RQ.entry[index].event_cycle < current_core_cycle[packet->cpu])
			RQ.entry[index].event_cycle = current_core_cycle[packet->cpu] + LATENCY;
		else
			RQ.entry[index].event_cycle += LATENCY;

		RQ.occupancy++;
		RQ.tail++;
		if (RQ.tail >= RQ.SIZE)
			RQ.tail = 0;

		DP ( if (warmup_complete[RQ.entry[index].cpu]) {
				cout << "[" << NAME << "_RQ] " <<  __func__ << " instr_id: " << RQ.entry[index].instr_id << " address: " << hex << RQ.entry[index].address;
				cout << " full_addr: " << RQ.entry[index].full_addr << dec;
				cout << " type: " << +RQ.entry[index].type << " head: " << RQ.head << " tail: " << RQ.tail << " occupancy: " << RQ.occupancy;
				cout << " event: " << RQ.entry[index].event_cycle << " current: " << current_core_cycle[RQ.entry[index].cpu] << endl; });

		if (packet->address == 0)
			assert(0);

		RQ.TO_CACHE++;
		RQ.ACCESS++;

		return -1;
	}

	int CACHE::add_wq(PACKET *packet)
	{
		// check for duplicates in the write queue
		int index = WQ.check_queue(packet);
		if (index != -1) {

			WQ.MERGED++;
			WQ.ACCESS++;

			return index; // merged index
		}

		// sanity check
		if (WQ.occupancy >= WQ.SIZE)
			assert(0);

		// if there is no duplicate, add it to the write queue
		index = WQ.tail;
		if (WQ.entry[index].address != 0) {
			cerr << "[" << NAME << "_ERROR] " << __func__ << " is not empty index: " << index;
			cerr << " address: " << hex << WQ.entry[index].address;
			cerr << " full_addr: " << WQ.entry[index].full_addr << dec << endl;
			assert(0);
		}

		WQ.entry[index] = *packet;

		// ADD LATENCY
		if (WQ.entry[index].event_cycle < current_core_cycle[packet->cpu])
			WQ.entry[index].event_cycle = current_core_cycle[packet->cpu] + LATENCY;
		else
			WQ.entry[index].event_cycle += LATENCY;

		WQ.occupancy++;
		WQ.tail++;
		if (WQ.tail >= WQ.SIZE)
			WQ.tail = 0;

		DP (if (warmup_complete[WQ.entry[index].cpu]) {
				cout << "[" << NAME << "_WQ] " <<  __func__ << " instr_id: " << WQ.entry[index].instr_id << " address: " << hex << WQ.entry[index].address;
				cout << " full_addr: " << WQ.entry[index].full_addr << dec;
				cout << " head: " << WQ.head << " tail: " << WQ.tail << " occupancy: " << WQ.occupancy;
				cout << " data: " << hex << WQ.entry[index].data << dec;
				cout << " event: " << WQ.entry[index].event_cycle << " current: " << current_core_cycle[WQ.entry[index].cpu] << endl; });

		WQ.TO_CACHE++;
		WQ.ACCESS++;

		return -1;
	}




	int CACHE::prefetch_page(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int fill_level, int pq_id, int free, int update_free, int free_distance, uint64_t id, int type, int iflag, int lad, int confidence, int irip)
	{
		int index, debug = 0, flag = 0, fctb_search = -10;
		uint64_t temp = va_to_pa_prefetch(cpu, base_addr, pf_addr), foo;

		if(!free)
			fctb_search = search_fctb(pf_addr);

		if(pq_id == 0){
			pf_requested++;
			pf_total_pq++;
		}

		if(pq_id == 0){
			index = PQ.check_queue_vpn(pf_addr).first;
		}
		else{
			cout << "I am using only one PQ" << endl;
		}

		if(temp){
			if(pq_id != 2){
				if(index != -1){
					if(debug)
						cout << "Duplicate in the Prefetch Queue: " << pf_addr << endl;
					return 0;
				}
			}

			PACKET pf_packet;
			pf_packet.data_pa = temp;

			if(pq_id == 0){
				if (PQ.occupancy == PQ.SIZE){
					uint64_t removed_vpn = PQ.remove_queue_lru();
				}
			}
			else{
				cout << "I am using only one PQ" << endl;
			}

			pf_packet.fill_level = fill_level;
			pf_packet.cpu = cpu;
			pf_packet.data = temp;
			pf_packet.address = pf_addr;
			pf_packet.full_addr = base_addr;
			pf_packet.ip = ip;
			pf_packet.type = 0xdeadbeef;
			pf_packet.event_cycle = current_core_cycle[cpu];
			pf_packet.free_bit = free;
			pf_packet.free_distance = free_distance;
			pf_packet.lad = lad;
			pf_packet.conf = confidence;
			pf_packet.irip = irip;

			if(fctb_search == -10){
				fctb_misses++;
				pf_packet.event_cycle = current_core_cycle[cpu];
				pf_packet.free_bit = free;
			}
			else{
				fctb_hits++;
				pf_packet.event_cycle = fctb[fctb_search][2];
				pf_packet.free_bit = 1;
			}

			if(free){
				if(lad == 0)
					pf_free++;
				pf_packet.stall_cycles = 100;
			}
			else{
				if(fctb_search != -10){
					pf_free++;
					pf_packet.stall_cycles = fctb[fctb_search][3];
				}
				else{
					pf_real++;
					int stall_cycles = mmu_cache_prefetch_search(cpu, pf_addr, 0, id, ip, type, iflag);
					pf_packet.stall_cycles = stall_cycles;

					int victim_entry = fctb_replacement_policy();
					fctb[victim_entry][0] = pf_addr;
					fctb[victim_entry][1] = (pf_addr & 0x07);
					fctb[victim_entry][2] = current_core_cycle[cpu];
					fctb[victim_entry][3] = stall_cycles;
				}
			}

			if(P2TLB){
				if(IS_STLB){
					uint32_t set = get_set(pf_addr);
					uint32_t way = get_way(pf_addr, set);

					way = find_victim(cpu, 0xdeadbeef, set, block[set], 0xdeadbeef, 0xdeadbeef, 0xdeadbeef);

					update_replacement_state(cpu, set, way, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, WRITEBACK, 1);

					block[set][way].valid = 1;

					block[set][way].dirty = 0;
					block[set][way].prefetch = 1;
					block[set][way].used = 0;

					block[set][way].tag = pf_addr;
					block[set][way].address = pf_addr;
					block[set][way].full_addr = 0xdeadbeef; 
					block[set][way].data = temp;
					block[set][way].cpu = cpu;
					block[set][way].stalls = pf_packet.event_cycle + pf_packet.stall_cycles;
				}

				return 1;
			}


			if(pq_id == 0)
				add_pq(&pf_packet);
			else
				cout << "I am using only one PQ" << endl;

			if(lad == 1)
				issued_prefetches_lad++;

			pf_issued++;
			return 1;
		}
		else{
			if(pq_id == 0)
				pf_swap++;
			return 0;
		}
	}

	int CACHE::prefetch_line(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int pf_fill_level, uint32_t prefetch_metadata)
	{
		pf_requested++;

		if (PQ.occupancy < PQ.SIZE) {
			if ((base_addr>>LOG2_PAGE_SIZE) == (pf_addr>>LOG2_PAGE_SIZE)) {

				PACKET pf_packet;
				pf_packet.fill_level = pf_fill_level;
				pf_packet.pf_origin_level = fill_level;
				if(pf_fill_level == FILL_L1)
				{
					pf_packet.fill_l1d = 1;
				}
				pf_packet.pf_metadata = prefetch_metadata;
				pf_packet.cpu = cpu;
				//pf_packet.data_index = LQ.entry[lq_index].data_index;
				//pf_packet.lq_index = lq_index;
				pf_packet.address = pf_addr >> LOG2_BLOCK_SIZE;
				pf_packet.full_addr = pf_addr;
				//pf_packet.instr_id = LQ.entry[lq_index].instr_id;
				//pf_packet.rob_index = LQ.entry[lq_index].rob_index;
				pf_packet.ip = ip;
				pf_packet.type = PREFETCH;
				pf_packet.event_cycle = current_core_cycle[cpu];

				// give a dummy 0 as the IP of a prefetch
				add_pq(&pf_packet);

				pf_issued++;

				return 1;
			}
		}

		return 0;
	}

	int CACHE::kpc_prefetch_line(uint64_t base_addr, uint64_t pf_addr, int pf_fill_level, int delta, int depth, int signature, int confidence, uint32_t prefetch_metadata)
	{
		if (PQ.occupancy < PQ.SIZE) {
			if ((base_addr>>LOG2_PAGE_SIZE) == (pf_addr>>LOG2_PAGE_SIZE)) {

				PACKET pf_packet;
				pf_packet.fill_level = pf_fill_level;
				pf_packet.pf_origin_level = fill_level;
				if(pf_fill_level == FILL_L1)
				{
					pf_packet.fill_l1d = 1;
				}
				pf_packet.pf_metadata = prefetch_metadata;
				pf_packet.cpu = cpu;
				//pf_packet.data_index = LQ.entry[lq_index].data_index;
				//pf_packet.lq_index = lq_index;
				pf_packet.address = pf_addr >> LOG2_BLOCK_SIZE;
				pf_packet.full_addr = pf_addr;
				//pf_packet.instr_id = LQ.entry[lq_index].instr_id;
				//pf_packet.rob_index = LQ.entry[lq_index].rob_index;
				pf_packet.ip = 0;
				pf_packet.type = PREFETCH;
				pf_packet.delta = delta;
				pf_packet.depth = depth;
				pf_packet.signature = signature;
				pf_packet.confidence = confidence;
				pf_packet.event_cycle = current_core_cycle[cpu];

				// give a dummy 0 as the IP of a prefetch
				add_pq(&pf_packet);

				pf_issued++;

				return 1;
			}
		}

		return 0;
	}

	int CACHE::add_pq(PACKET *packet)
	{
		// check for the latest wirtebacks in the write queue
		int wq_index = WQ.check_queue(packet);
		if (wq_index != -1) {

			// check fill level
			if (packet->fill_level < fill_level) {

				packet->data = WQ.entry[wq_index].data;

				if(fill_level == FILL_L2)
				{
					if(packet->fill_l1i)
					{
						upper_level_icache[packet->cpu]->return_data(packet);
					}
					if(packet->fill_l1d)
					{
						upper_level_dcache[packet->cpu]->return_data(packet);
					}
				}
				else
				{
					if (packet->instruction)
						upper_level_icache[packet->cpu]->return_data(packet);
					if (packet->is_data)
						upper_level_dcache[packet->cpu]->return_data(packet);
				}
			}

			HIT[packet->type]++;
			ACCESS[packet->type]++;

			WQ.FORWARD++;
			PQ.ACCESS++;

			return -1;
		}

		// check for duplicates in the PQ

		int index;
		if (cache_type == IS_STLB)
			index = PQ.check_queue_vpn(packet->address).first;
		else
			index = PQ.check_queue(packet);

		if (index != -1) {
			if (packet->fill_level < PQ.entry[index].fill_level)
			{
				PQ.entry[index].fill_level = packet->fill_level;
			}
			if((packet->instruction == 1) && (PQ.entry[index].instruction != 1))
			{
				PQ.entry[index].instruction = 1;
			}
			if((packet->is_data == 1) && (PQ.entry[index].is_data != 1))
			{
				PQ.entry[index].is_data = 1;
			}
			if((packet->fill_l1i) && (PQ.entry[index].fill_l1i != 1))
			{
				PQ.entry[index].fill_l1i = 1;
			}
			if((packet->fill_l1d) && (PQ.entry[index].fill_l1d != 1))
			{
				PQ.entry[index].fill_l1d = 1;
			}

			PQ.MERGED++;
			PQ.ACCESS++;

			return index; // merged index
		}

		// check occupancy
		if (PQ.occupancy == PQ_SIZE) {
			PQ.FULL++;

			DP ( if (warmup_complete[packet->cpu]) {
					cout << "[" << NAME << "] cannot process add_pq since it is full" << endl; });
			return -2; // cannot handle this request
		}

		/* if there is no duplicate, add it to PQ */
		if(cache_type != IS_STLB)
			index = PQ.tail;
		else{
			while(1){
				int test = PQ.update_tail();
				if(test == 64){
					PQ.remove_queue_lru();
					test = PQ.update_tail();
				}
				else
					break;
			}

			index = PQ.tail;
		}

		if(cache_type == IS_STLB){
			//cout << hex << PQ.entry[index].address << endl;
			//cout << dec << PQ.occupancy << endl;
			//cout << "--------------------------------" << endl;
			//PQ.print_pq();
		}
#ifdef SANITY_CHECK

		if (PQ.entry[index].address != 0) {
			cerr << "[" << NAME << "_ERROR] " << __func__ << " is not empty index: " << index;
			cerr << " address: " << hex << PQ.entry[index].address;
			cerr << " full_addr: " << PQ.entry[index].full_addr << dec << endl;
			assert(0);
		}
#endif

		PQ.entry[index] = *packet;

		// ADD LATENCY
		if (PQ.entry[index].event_cycle < current_core_cycle[packet->cpu])
			PQ.entry[index].event_cycle = current_core_cycle[packet->cpu] + LATENCY;
		else
			PQ.entry[index].event_cycle += LATENCY;

		PQ.occupancy++;
		PQ.tail++;
		if (PQ.tail >= PQ.SIZE)
			PQ.tail = 0;

		DP ( if (warmup_complete[PQ.entry[index].cpu]) {
				cout << "[" << NAME << "_PQ] " <<  __func__ << " instr_id: " << PQ.entry[index].instr_id << " address: " << hex << PQ.entry[index].address;
				cout << " full_addr: " << PQ.entry[index].full_addr << dec;
				cout << " type: " << +PQ.entry[index].type << " head: " << PQ.head << " tail: " << PQ.tail << " occupancy: " << PQ.occupancy;
				cout << " event: " << PQ.entry[index].event_cycle << " current: " << current_core_cycle[PQ.entry[index].cpu] << endl; });

		if (packet->address == 0)
			assert(0);

		PQ.TO_CACHE++;
		PQ.ACCESS++;

		return -1;
	}

	void CACHE::return_data(PACKET *packet)
	{
		// check MSHR information
		int mshr_index = check_mshr(packet);

		// sanity check
		if (mshr_index == -1) {
			cerr << "[" << NAME << "_MSHR] " << __func__ << " instr_id: " << packet->instr_id << " cannot find a matching entry!";
			cerr << " full_addr: " << hex << packet->full_addr;
			cerr << " address: " << packet->address << dec;
			cerr << " event: " << packet->event_cycle << " current: " << current_core_cycle[packet->cpu] << endl;
			assert(0);
		}

		// MSHR holds the most updated information about this request
		// no need to do memcpy
		MSHR.num_returned++;
		MSHR.entry[mshr_index].returned = COMPLETED;
		MSHR.entry[mshr_index].data = packet->data;
		MSHR.entry[mshr_index].pf_metadata = packet->pf_metadata;

		// ADD LATENCY
		if (MSHR.entry[mshr_index].event_cycle < current_core_cycle[packet->cpu])
			MSHR.entry[mshr_index].event_cycle = current_core_cycle[packet->cpu] + LATENCY;
		else
			MSHR.entry[mshr_index].event_cycle += LATENCY;

		update_fill_cycle();

		DP (if (warmup_complete[packet->cpu]) {
				cout << "[" << NAME << "_MSHR] " <<  __func__ << " instr_id: " << MSHR.entry[mshr_index].instr_id;
				cout << " address: " << hex << MSHR.entry[mshr_index].address << " full_addr: " << MSHR.entry[mshr_index].full_addr;
				cout << " data: " << MSHR.entry[mshr_index].data << dec << " num_returned: " << MSHR.num_returned;
				cout << " index: " << mshr_index << " occupancy: " << MSHR.occupancy;
				cout << " event: " << MSHR.entry[mshr_index].event_cycle << " current: " << current_core_cycle[packet->cpu] << " next: " << MSHR.next_fill_cycle << endl; });
	}

	void CACHE::update_fill_cycle()
	{
		// update next_fill_cycle
		uint64_t min_cycle = UINT64_MAX;
		uint32_t min_index = MSHR.SIZE;
		for (uint32_t i=0; i<MSHR.SIZE; i++) {
			if ((MSHR.entry[i].returned == COMPLETED) && (MSHR.entry[i].event_cycle < min_cycle)) {
				min_cycle = MSHR.entry[i].event_cycle;
				min_index = i;
			}

			DP (if (warmup_complete[MSHR.entry[i].cpu]) {
					cout << "[" << NAME << "_MSHR] " <<  __func__ << " checking instr_id: " << MSHR.entry[i].instr_id;
					cout << " address: " << hex << MSHR.entry[i].address << " full_addr: " << MSHR.entry[i].full_addr;
					cout << " data: " << MSHR.entry[i].data << dec << " returned: " << +MSHR.entry[i].returned << " fill_level: " << MSHR.entry[i].fill_level;
					cout << " index: " << i << " occupancy: " << MSHR.occupancy;
					cout << " event: " << MSHR.entry[i].event_cycle << " current: " << current_core_cycle[MSHR.entry[i].cpu] << " next: " << MSHR.next_fill_cycle << endl; });
		}

		MSHR.next_fill_cycle = min_cycle;
		MSHR.next_fill_index = min_index;
		if (min_index < MSHR.SIZE) {

			DP (if (warmup_complete[MSHR.entry[min_index].cpu]) {
					cout << "[" << NAME << "_MSHR] " <<  __func__ << " instr_id: " << MSHR.entry[min_index].instr_id;
					cout << " address: " << hex << MSHR.entry[min_index].address << " full_addr: " << MSHR.entry[min_index].full_addr;
					cout << " data: " << MSHR.entry[min_index].data << dec << " num_returned: " << MSHR.num_returned;
					cout << " event: " << MSHR.entry[min_index].event_cycle << " current: " << current_core_cycle[MSHR.entry[min_index].cpu] << " next: " << MSHR.next_fill_cycle << endl; });
		}
	}

	int CACHE::check_mshr(PACKET *packet)
	{
		// search mshr
		//bool instruction_and_data_collision = false;

		for (uint32_t index=0; index<MSHR_SIZE; index++)
		{
			if (MSHR.entry[index].address == packet->address)
			{
				//if(MSHR.entry[index].instruction != packet->instruction)
				//  {
				//    instruction_and_data_collision = true;
				//  }
				//else
				//  {
				DP ( if (warmup_complete[packet->cpu]) {
						cout << "[" << NAME << "_MSHR] " << __func__ << " same entry instr_id: " << packet->instr_id << " prior_id: " << MSHR.entry[index].instr_id;
						cout << " address: " << hex << packet->address;
						cout << " full_addr: " << packet->full_addr << dec << endl; });

				return index;
				//  }
			}
		}

		//if(instruction_and_data_collision) // remove instruction-and-data collision safeguard
		//  {
		//return -2;
		//  }

		DP ( if (warmup_complete[packet->cpu]) {
				cout << "[" << NAME << "_MSHR] " << __func__ << " new address: " << hex << packet->address;
				cout << " full_addr: " << packet->full_addr << dec << endl; });

		DP ( if (warmup_complete[packet->cpu] && (MSHR.occupancy == MSHR_SIZE)) { 
				cout << "[" << NAME << "_MSHR] " << __func__ << " mshr is full";
				cout << " instr_id: " << packet->instr_id << " mshr occupancy: " << MSHR.occupancy;
				cout << " address: " << hex << packet->address;
				cout << " full_addr: " << packet->full_addr << dec;
				cout << " cycle: " << current_core_cycle[packet->cpu] << endl; });

		return -1;
	}

	void CACHE::add_mshr(PACKET *packet)
	{
		uint32_t index = 0;

		packet->cycle_enqueued = current_core_cycle[packet->cpu];

		// search mshr
		for (index=0; index<MSHR_SIZE; index++) {
			if (MSHR.entry[index].address == 0) {

				MSHR.entry[index] = *packet;
				MSHR.entry[index].returned = INFLIGHT;
				MSHR.occupancy++;

				DP ( if (warmup_complete[packet->cpu]) {
						cout << "[" << NAME << "_MSHR] " << __func__ << " instr_id: " << packet->instr_id;
						cout << " address: " << hex << packet->address << " full_addr: " << packet->full_addr << dec;
						cout << " index: " << index << " occupancy: " << MSHR.occupancy << endl; });

				break;
			}
		}
	}

	uint32_t CACHE::get_occupancy(uint8_t queue_type, uint64_t address)
	{
		if (queue_type == 0)
			return MSHR.occupancy;
		else if (queue_type == 1)
			return RQ.occupancy;
		else if (queue_type == 2)
			return WQ.occupancy;
		else if (queue_type == 3)
			return PQ.occupancy;

		return 0;
	}

	uint32_t CACHE::get_size(uint8_t queue_type, uint64_t address)
	{
		if (queue_type == 0)
			return MSHR.SIZE;
		else if (queue_type == 1)
			return RQ.SIZE;
		else if (queue_type == 2)
			return WQ.SIZE;
		else if (queue_type == 3)
			return PQ.SIZE;

		return 0;
	}

	void CACHE::increment_WQ_FULL(uint64_t address)
	{
		WQ.FULL++;
	}
