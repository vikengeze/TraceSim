#ifndef CACHE_H
#define CACHE_H

#include "memory_class.h"
#include <map>
#include <iterator>

// PAGE
extern uint32_t PAGE_TABLE_LATENCY, SWAP_LATENCY;

#define P2TLB 1

// Free Prefetching

// flag for demand page walks
#define ENABLE_FP 0

// flag for prefetch page walks
#define ENABLE_PREF_FP 1

// Lookahead Depth (0:disable)
#define LA_DEPTH 0

// Replacement Policy for Markov's prediction table --> [ 0:LRU, 1:LFU, 2:RANDOM, 3:... ]
#define RP_MP 1

// Number of prediction table entries you randomly select from for eviction in LFU replacement policy for the Markov instruction TLB prefetcher
#define LLIMIT 5

// Number of bits for the confidence counters of the RP_SUC_MP replacement policy
#define CNF_BITS 3

// Number of successors of Markov I-TLB Prefetcher
#define SUCCESSORS 2

// Replacement Policy for the successors of Markov's prediction table --> [ 0:FIFO, 1:RANDOM, 2:CUSTOM, 3:... ]
#define RP_SUC_MP 2

// Number of STLB instruction misses to reset frequency field of Markov Prefetcher used for the LFU policy
#define RESET_FREQ 5000

#define PML4_SET 2
#define PDP_SET 4
#define PD_SET 16

#define PML4_WAY 1
#define PDP_WAY 1
#define PD_WAY 2

// CloudSuite
#define CLOUD_SUITE 0

#define FCTB_SIZE 4

// CACHE TYPE
#define IS_ITLB 0
#define IS_DTLB 1
#define IS_STLB 2
#define IS_L1I  3
#define IS_L1D  4
#define IS_L2C  5
#define IS_LLC  6

// INSTRUCTION TLB
#define ITLB_SET 16
#define ITLB_WAY 8
#define ITLB_RQ_SIZE 16
#define ITLB_WQ_SIZE 16
#define ITLB_PQ_SIZE 0
#define ITLB_MSHR_SIZE 4
#define ITLB_LATENCY 1

// DATA TLB
#define DTLB_SET 16
#define DTLB_WAY 4
#define DTLB_RQ_SIZE 16
#define DTLB_WQ_SIZE 16
#define DTLB_PQ_SIZE 0
#define DTLB_MSHR_SIZE 4
#define DTLB_LATENCY 1

// SECOND LEVEL TLB
#define STLB_SET 256
#define STLB_WAY 6
#define STLB_RQ_SIZE 32
#define STLB_WQ_SIZE 32
#define STLB_PQ_SIZE 64
#define STLB_MSHR_SIZE 4
#define STLB_LATENCY 8

// L1 INSTRUCTION CACHE
#define L1I_SET 64
#define L1I_WAY 8
#define L1I_RQ_SIZE 64
#define L1I_WQ_SIZE 64 
#define L1I_PQ_SIZE 32
#define L1I_MSHR_SIZE 8
#define L1I_LATENCY 4

// L1 DATA CACHE
#define L1D_SET 64
#define L1D_WAY 8
#define L1D_RQ_SIZE 64
#define L1D_WQ_SIZE 64 
#define L1D_PQ_SIZE 8
#define L1D_MSHR_SIZE 16 //8 //16
#define L1D_LATENCY 4

// L2 CACHE
#define L2C_SET 1024 //512
#define L2C_WAY 8
#define L2C_RQ_SIZE 32
#define L2C_WQ_SIZE 32
#define L2C_PQ_SIZE 32
#define L2C_MSHR_SIZE 32 //16 //32
#define L2C_LATENCY 8  // 4/5 (L1I or L1D) + 10 = 14/15 cycles

// LAST LEVEL CACHE
#define LLC_SET NUM_CPUS*2048
#define LLC_WAY 16
#define LLC_RQ_SIZE 48//NUM_CPUS*L2C_MSHR_SIZE //48
#define LLC_WQ_SIZE 48//NUM_CPUS*L2C_MSHR_SIZE //48
#define LLC_PQ_SIZE NUM_CPUS*64//1
#define LLC_MSHR_SIZE NUM_CPUS*64//32
#define LLC_LATENCY 10 //10  // 4/5 (L1I or L1D) + 10 + 20 = 34/35 cycles

class CACHE : public MEMORY {
	public:
		uint32_t cpu;
		const string NAME;
		const uint32_t NUM_SET, NUM_WAY, NUM_LINE, WQ_SIZE, RQ_SIZE, PQ_SIZE, MSHR_SIZE;
		uint32_t LATENCY;
		BLOCK **block;
		int fill_level;
		uint32_t MAX_READ, MAX_FILL;
		uint32_t reads_available_this_cycle;
		uint8_t cache_type;

		uint64_t fctb[FCTB_SIZE][4];

		uint64_t timer;

		int selector;

		map<uint64_t, uint64_t> code_footprint;
		map<uint64_t, uint64_t>::iterator it;

		// prefetch stats
		uint64_t pf_requested,
			 pf_issued,
			 pf_useful,
			 pf_useless,
			 pf_fill,
			 pf_hits_pq,
			 pf_misses_pq,
			 pf_swap,
			 pf_dupli,
			 pf_free,
			 pf_real,
			 previous_iva,
			 previous_ip,
			 fctb_hits,
			 fctb_misses, 
			 hit_prefetches_lad,
			 issued_prefetches_lad,
			 pf_total_pq,
			 morrigan_filter_hits,
			 irip_hits,
			 sdp_hits,
			 bpbp[5]; // 0: # instruction prefetches 1: portion of instruction prefetches that are in the same page 2: portion of beyond page boundaries instruction prefetches 3: beyond page boundaries prefetches that hit in the TLB hierarchy 4: beyond page boundaries prefetches that miss in the TLB

		// queues
		PACKET_QUEUE WQ{NAME + "_WQ", WQ_SIZE}, // write queue
			     RQ{NAME + "_RQ", RQ_SIZE}, // read queue
			     PQ{NAME + "_PQ", PQ_SIZE}, // prefetch queue
			     MSHR{NAME + "_MSHR", MSHR_SIZE}, // MSHR
			     PROCESSED{NAME + "_PROCESSED", ROB_SIZE}; // processed queue

		uint64_t sim_access[NUM_CPUS][NUM_TYPES],
		sim_hit[NUM_CPUS][NUM_TYPES],
		sim_miss[NUM_CPUS][NUM_TYPES],
		roi_access[NUM_CPUS][NUM_TYPES],
		roi_hit[NUM_CPUS][NUM_TYPES],
		roi_miss[NUM_CPUS][NUM_TYPES];

		uint64_t free_distance_table[14];

		uint64_t pml4[PML4_SET][PML4_WAY], pdp[PDP_SET][PDP_WAY], pd[PD_SET][PD_WAY];
		uint64_t pml4_lru[PML4_SET][PML4_WAY], pdp_lru[PDP_SET][PDP_WAY], pd_lru[PD_SET][PD_WAY];
		uint64_t mmu_cache_demand_hits[4], mmu_cache_prefetch_hits[4];
		uint64_t mmu_timer;

		uint64_t pagetable_mr_hit_ratio[4][4];
		uint64_t rfhits[2];
		uint64_t free_hits[14];

		uint64_t decay_timer;
		uint64_t decay_conf_timer;

		uint64_t stlb_misses[2]; // 0 is data, 1 is instruction

		uint64_t total_miss_latency;

		// constructor
		CACHE(string v1, uint32_t v2, int v3, uint32_t v4, uint32_t v5, uint32_t v6, uint32_t v7, uint32_t v8) 
			: NAME(v1), NUM_SET(v2), NUM_WAY(v3), NUM_LINE(v4), WQ_SIZE(v5), RQ_SIZE(v6), PQ_SIZE(v7), MSHR_SIZE(v8) {

				LATENCY = 0;

				// cache block
				block = new BLOCK* [NUM_SET];
				for (uint32_t i=0; i<NUM_SET; i++) {
					block[i] = new BLOCK[NUM_WAY]; 

					for (uint32_t j=0; j<NUM_WAY; j++) {
						block[i][j].lru = j;
					}
				}

				// fctb initialization loop
				for (uint32_t j=0; j<FCTB_SIZE; j++) {
					fctb[j][0] = 0;
					fctb[j][1] = 0;
					fctb[j][2] = 0;
					fctb[j][3] = 0;
				}

				for (int i=0; i<PML4_SET; i++){
					for(int j=0; j<PML4_WAY; j++){
						pml4[i][j] = 0;
						pml4_lru[i][j] = 0;
					}
				}

				for (int i=0; i<PDP_SET; i++){
					for(int j=0; j<PDP_WAY; j++){
						pdp[i][j] = 0;
						pdp_lru[i][j] = 0;
					}
				}

				for (int i=0; i<PD_SET; i++){
					for(int j=0; j<PD_WAY; j++){
						pd[i][j] = 0;
						pd_lru[i][j] = 0;
					}
				}

				mmu_timer = 0;

				for(int i=0; i<4; i++){
					mmu_cache_demand_hits[i] = 0;
					mmu_cache_prefetch_hits[i] = 0;
				}

				for(int i=0; i<4; i++){
					for(int j=0; j<4; j++){
						pagetable_mr_hit_ratio[i][j] = 0;
					}
				}

				for(int i=0; i<14; i++){
					free_hits[i] = 0;
				}

				rfhits[0] = 0;
				rfhits[1] = 0;
				decay_timer = 0;
				timer = 0;

				pf_hits_pq = 0;
				pf_misses_pq = 0;
				pf_swap = 0;
				pf_dupli = 0;
				pf_free = 0;
				pf_real = 0;
				previous_iva = 0;
				previous_ip = 0;
				stlb_misses[0] = 0;
				stlb_misses[1] = 0;
				fctb_hits = 0;
				fctb_misses = 0;
				issued_prefetches_lad = 0;
				hit_prefetches_lad = 0;
				pf_total_pq = 0;

				selector = 0;

				irip_hits = 0;
				sdp_hits = 0;

				bpbp[0] = 0;
				bpbp[1] = 0;
				bpbp[2] = 0;
				bpbp[3] = 0;
				bpbp[4] = 0;

				morrigan_filter_hits = 0;

				for (uint32_t i=0; i<NUM_CPUS; i++) {
					upper_level_icache[i] = NULL;
					upper_level_dcache[i] = NULL;

					for (uint32_t j=0; j<NUM_TYPES; j++) {
						sim_access[i][j] = 0;
						sim_hit[i][j] = 0;
						sim_miss[i][j] = 0;
						roi_access[i][j] = 0;
						roi_hit[i][j] = 0;
						roi_miss[i][j] = 0;
					}
				}

				total_miss_latency = 0;

				lower_level = NULL;
				extra_interface = NULL;
				fill_level = -1;
				MAX_READ = 1;
				MAX_FILL = 1;

				pf_requested = 0;
				pf_issued = 0;
				pf_useful = 0;
				pf_useless = 0;
				pf_fill = 0;
			};

		// destructor
		~CACHE() {
			for (uint32_t i=0; i<NUM_SET; i++)
				delete[] block[i];
			delete[] block;
		};

		// functions
		pair<int, int> check_hit_stlb_pq(uint64_t vpn);

		int  add_rq(PACKET *packet),
		     add_wq(PACKET *packet),
		     add_pq(PACKET *packet);

		void return_data(PACKET *packet),
		     operate(),
		     increment_WQ_FULL(uint64_t address);

		uint32_t get_occupancy(uint8_t queue_type, uint64_t address),
			 get_size(uint8_t queue_type, uint64_t address);

		int  check_hit(PACKET *packet),
		     invalidate_entry(uint64_t inval_addr),
		     check_mshr(PACKET *packet),
		     prefetch_line(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int prefetch_fill_level, uint32_t prefetch_metadata),
		     kpc_prefetch_line(uint64_t base_addr, uint64_t pf_addr, int prefetch_fill_level, int delta, int depth, int signature, int confidence, uint32_t prefetch_metadata),
		     prefetch_page(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int fill_level, int pq_id, int free, int update_free, int free_distance, uint64_t id, int type, int iflag, int lad, int confidence, int irip);

		void handle_fill(),
		     handle_writeback(),
		     handle_read(),
		     handle_prefetch();

		void add_mshr(PACKET *packet),
		     update_fill_cycle(),
		     llc_initialize_replacement(),
		     update_replacement_state(uint32_t cpu, uint32_t set, uint32_t way, uint64_t full_addr, uint64_t ip, uint64_t victim_addr, uint32_t type, uint8_t hit),
		     llc_update_replacement_state(uint32_t cpu, uint32_t set, uint32_t way, uint64_t full_addr, uint64_t ip, uint64_t victim_addr, uint32_t type, uint8_t hit),
		     lru_update(uint32_t set, uint32_t way),
		     fill_cache(uint32_t set, uint32_t way, PACKET *packet),
		     replacement_final_stats(),
		     llc_replacement_final_stats(),
		     //prefetcher_initialize(),
		     l1d_prefetcher_initialize(),
		     l2c_prefetcher_initialize(),
		     llc_prefetcher_initialize(),
		     stlb_prefetcher_initialize(),
		     prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type),
		     l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type),
		     prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr),
		     stlb_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr),
		     l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in),
		     //prefetcher_final_stats(),
		     l1d_prefetcher_final_stats(),
		     l2c_prefetcher_final_stats(),
		     llc_prefetcher_final_stats(),
		     stlb_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, int answer, int warmup, int * free_indexes, uint64_t instr_id, int iflag, uint64_t curr_cycle),
		     stlb_prefetcher_final_stats(uint64_t prefetches, uint64_t hits, uint64_t misses, uint64_t swap, uint64_t dupli, uint64_t free, uint64_t real, uint64_t *mmu_cache_demand_hits, uint64_t * mmu_cache_prefetch_hits, uint64_t * rfhits, uint64_t * free_hits, uint64_t mr[4][4], uint64_t stlb_misses[2]);

		int  search_pml4(uint64_t address), search_pdp(uint64_t address), search_pd(uint64_t address);
		void lru_pml4(uint64_t timer, uint64_t address), lru_pdp(uint64_t timer, uint64_t address), lru_pd(uint64_t timer, uint64_t address);
		void print_pml4(), print_pdp(), print_pd();
		void free_prefetching(uint64_t ip, uint64_t addr, int cache_line_position_n, uint64_t pf_addr, int * free_indexes, uint64_t instr_id, int type, int iflag);

		// lookahead prefetching for Markov I-TLB Prefetcher
		void lookahead_prefetching(uint64_t ip, uint64_t addr, uint64_t pf_addr, uint64_t instr_id, int type, int iflag, int depth);

		// FCTB functions
		int fctb_replacement_policy(), search_fctb(uint64_t current_vpn);
		void print_fctb(), refresh_fctb(uint64_t current_cycle);

		void (*l1i_prefetcher_cache_operate)(uint32_t, uint64_t, uint8_t, uint8_t);
		void (*l1i_prefetcher_cache_fill)(uint32_t, uint64_t, uint32_t, uint32_t, uint8_t, uint64_t);

		uint32_t l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in),
			 llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in),
			 l2c_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in),
			 llc_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in);

		uint32_t get_set(uint64_t address),
			 get_way(uint64_t address, uint32_t set),
			 find_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type),
			 llc_find_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type),
			 lru_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type);

		int * sorted_free_distances();

		void issue_prefetches(int offset, uint64_t current_vpn, uint64_t ip, uint64_t instr_id, int iflag, int * free_indexes, uint8_t type, int answer);
};

#endif
