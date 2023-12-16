//============================================================================
// 
//  SAM Coupe replica for MiST
//  Copyright (C) 2016-2018 Sorgelig
//
//  This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the Free
//  Software Foundation; either version 2 of the License, or (at your option)
//  any later version.
//
//  This program is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//  more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
//============================================================================

`default_nettype none

module SamCoupe
(
	input         CLOCK_27,
`ifdef USE_CLOCK_50
	input         CLOCK_50,
`endif
	output        LED,
	output [VGA_BITS-1:0] VGA_R,
	output [VGA_BITS-1:0] VGA_G,
	output [VGA_BITS-1:0] VGA_B,
	output        VGA_HS,
	output        VGA_VS,

`ifdef USE_HDMI
	output        HDMI_RST,
	output  [7:0] HDMI_R,
	output  [7:0] HDMI_G,
	output  [7:0] HDMI_B,
	output        HDMI_HS,
	output        HDMI_VS,
	output        HDMI_PCLK,
	output        HDMI_DE,
	inout         HDMI_SDA,
	inout         HDMI_SCL,
	input         HDMI_INT,
`endif

	input         SPI_SCK,
	inout         SPI_DO,
	input         SPI_DI,
	input         SPI_SS2,    // data_io
	input         SPI_SS3,    // OSD
	input         CONF_DATA0, // SPI_SS for user_io

`ifdef USE_QSPI
	input         QSCK,
	input         QCSn,
	inout   [3:0] QDAT,
`endif
`ifndef NO_DIRECT_UPLOAD
	input         SPI_SS4,
`endif

	output [12:0] SDRAM_A,
	inout  [15:0] SDRAM_DQ,
	output        SDRAM_DQML,
	output        SDRAM_DQMH,
	output        SDRAM_nWE,
	output        SDRAM_nCAS,
	output        SDRAM_nRAS,
	output        SDRAM_nCS,
	output  [1:0] SDRAM_BA,
	output        SDRAM_CLK,
	output        SDRAM_CKE,

`ifdef DUAL_SDRAM
	output [12:0] SDRAM2_A,
	inout  [15:0] SDRAM2_DQ,
	output        SDRAM2_DQML,
	output        SDRAM2_DQMH,
	output        SDRAM2_nWE,
	output        SDRAM2_nCAS,
	output        SDRAM2_nRAS,
	output        SDRAM2_nCS,
	output  [1:0] SDRAM2_BA,
	output        SDRAM2_CLK,
	output        SDRAM2_CKE,
`endif

	output        AUDIO_L,
	output        AUDIO_R,
`ifdef I2S_AUDIO
	output        I2S_BCK,
	output        I2S_LRCK,
	output        I2S_DATA,
`endif
`ifdef SPDIF_AUDIO
	output        SPDIF,
`endif
`ifdef USE_AUDIO_IN
	input         AUDIO_IN,
`endif
	input         UART_RX,
	output        UART_TX
);

`ifdef NO_DIRECT_UPLOAD
localparam bit DIRECT_UPLOAD = 0;
wire SPI_SS4 = 1;
`else
localparam bit DIRECT_UPLOAD = 1;
`endif

`ifdef USE_QSPI
localparam bit QSPI = 1;
assign QDAT = 4'hZ;
`else
localparam bit QSPI = 0;
`endif

`ifdef VGA_8BIT
localparam VGA_BITS = 8;
`else
localparam VGA_BITS = 6;
`endif

`ifdef USE_HDMI
localparam bit HDMI = 1;
assign HDMI_RST = 1'b1;
`else
localparam bit HDMI = 0;
`endif

`ifdef BIG_OSD
localparam bit BIG_OSD = 1;
`define SEP "-;",
`else
localparam bit BIG_OSD = 0;
`define SEP
`endif

// remove this if the 2nd chip is actually used
`ifdef DUAL_SDRAM
assign SDRAM2_A = 13'hZZZZ;
assign SDRAM2_BA = 0;
assign SDRAM2_DQML = 0;
assign SDRAM2_DQMH = 0;
assign SDRAM2_CKE = 0;
assign SDRAM2_CLK = 0;
assign SDRAM2_nCS = 1;
assign SDRAM2_DQ = 16'hZZZZ;
assign SDRAM2_nCAS = 1;
assign SDRAM2_nRAS = 1;
assign SDRAM2_nWE = 1;
`endif

`include "build_id.v"

assign LED = ~ioctl_download;

localparam CONF_STR = 
{
	"SAMCOUPE;;",
	"S0,DSKMGTIMG,Drive 1;",
	"S1,DSKMGTIMG,Drive 2;",
	`SEP
	"O4,Drive Write,Prohibit,Allow;",
	"O12,Scandoubler Fx,None,CRT 25%,CRT 50%,CRT 75%;",
	"O8A,CPU Speed,Normal,6MHz,9.6MHz,12MHz,24MHz;",
	"OBC,ZX Mode Speed,Emulated,Full,Real;",
	"O5,External RAM,on,off;",
	"V,v1.56.",`BUILD_DATE
};

wire  [1:0] st_scanlines = status[2:1];

////////////////////   CLOCKS   ///////////////////
wire clk_sys, clk_mem;
wire locked;

assign SDRAM_CLK = clk_mem;
pll pll
(
`ifdef USE_CLOCK_50
	.inclk0(CLOCK_50),
`else
	.inclk0(CLOCK_27),
`endif
	.c0(clk_mem),
	.c2(clk_sys),
	.locked(locked)
);

reg  ce_1m;
reg  ce_psg;   //8MHz
reg  ce_6mp;
reg  ce_6mn;
reg  ce_24m;
reg  cpu_en;
reg  cpu_p;
reg  cpu_n;

wire ce_cpu_p = cpu_en & cpu_p;
wire ce_cpu_n = cpu_en & cpu_n;

reg [2:0] req_speed, cpu_speed = 0;
always_comb begin
	casex({turbo_boot, !video_mode && (status[12:11] == 2)})
		'b1X: req_speed = 5;
		'b01: req_speed = 0;
		'b00: req_speed = 3'd1 + status[10:8];
	default: req_speed = 1;
	endcase
end

reg turbo_boot = 0;
always @(posedge clk_sys) begin
	reg old_rd, skip;
	old_rd <= port_rd;

	if(reset) {skip,turbo_boot} <= 1;
	else if(~old_rd & port_rd & kbdr_sel) begin
		skip <= 1;
		if(skip) turbo_boot <= 0;
	end
end

wire [4:0] cpu_max[6] = '{13, 7, 7, 4, 3, 1};
wire [4:0] cpu_mid[6] = '{6,  4, 4, 2, 2, 0};

always @(posedge clk_sys) begin
	reg [2:0] counter = 0;
	reg [2:0] psg_div = 0;
	reg [4:0] cpu_div = 0;
	reg [5:0] meg_div = 0;
	reg       cnt_en  = 1;
	reg       skip;

	counter <=  counter + 1'd1;
	ce_24m  <= !counter[1];
	ce_6mp  <= !counter[2] & !counter[1:0];
	ce_6mn  <=  counter[2] & !counter[1:0];

	{cpu_p, cpu_n} <= 0;
	if(cnt_en & ~(ram_busy & (cpu_speed >= 3) & (cpu_div == (cpu_mid[cpu_speed])))) begin
		cpu_div <= cpu_div + 1'd1;
		if(cpu_div >= cpu_max[cpu_speed]) cpu_div <= 0;

		if(!cpu_div) cpu_en <= ~((cpu_speed == 1) & (ram_wait | io_wait));

		cpu_p <= (cpu_div == 0);
		cpu_n <= (cpu_div == cpu_mid[cpu_speed]);
	end

	if((cpu_speed != req_speed) & nWR & nRD & ce_cpu_p) begin
		cpu_div <= 1;
		cnt_en  <= 0;
		skip    <= 0;
	end

	if(~cnt_en & !counter) begin
		skip <= 1;
		if(skip) begin
			cpu_speed <= req_speed;
			cnt_en    <= 1;
		end
	end

	psg_div <= psg_div + 1'd1;
	if(psg_div == 5) psg_div <= 0;
	ce_psg  <= !psg_div;

	meg_div <= meg_div + 1'd1;
	if(meg_div == 47) meg_div <= 0;
	ce_1m  <= !meg_div;
end

// Contention model
wire ram_acc = ~nMREQ & nRFSH & ~rom0_sel & ~rom1_sel & ~ext_ram;
wire io_acc  = ~nIORQ & ~(nRD & nWR) & nM1 & &addr[7:3];
reg  ram_wait, io_wait;

always @(posedge clk_sys) begin
	reg old_ram, old_io, old_memcont, old_iocont;

	old_ram <= ram_acc;
	if(~old_ram & ram_acc & mem_contention) ram_wait <= 1;

	old_memcont <= mem_contention;
	if(~mem_contention & old_memcont) ram_wait <= 0;

	old_io  <= io_acc;
	if(~old_io & io_acc & io_contention) io_wait <= 1;

	old_iocont  <= io_contention;
	if(~io_contention & old_iocont) io_wait <= 0;
end


//////////////////   MIST ARM I/O   ///////////////////
wire        key_strobe;
wire        key_pressed;
wire        key_extended;
wire  [7:0] key_code;

wire  [8:0] mouse_x;
wire  [8:0] mouse_y;
wire  [7:0] mouse_flags;
wire        mouse_strobe;

wire [24:0] ps2_mouse = { mouse_strobe_level, mouse_y[7:0], mouse_x[7:0], mouse_flags };
reg         mouse_strobe_level;
always @(posedge clk_sys) if (mouse_strobe) mouse_strobe_level <= ~mouse_strobe_level;

wire  [7:0] joystick_0;
wire  [7:0] joystick_1;
wire  [1:0] buttons;
wire  [1:0] switches;
wire        scandoubler_disable;
wire        ypbpr;
wire        no_csync;
wire [31:0] status;

wire [31:0] sd_lba;
wire  [1:0] sd_rd;
wire  [1:0] sd_wr;
wire        sd_ack;
wire  [8:0] sd_buff_addr;
wire  [7:0] sd_buff_dout;
wire  [7:0] sd_buff_din;
wire        sd_buff_wr;
wire  [1:0] img_mounted;
wire [31:0] img_size;

`ifdef USE_HDMI
wire        i2c_start;
wire        i2c_read;
wire  [6:0] i2c_addr;
wire  [7:0] i2c_subaddr;
wire  [7:0] i2c_dout;
wire  [7:0] i2c_din;
wire        i2c_ack;
wire        i2c_end;
`endif

user_io #(.STRLEN($size(CONF_STR)>>3), .SD_IMAGES(2), .FEATURES(32'h0 | (BIG_OSD << 13) | (HDMI << 14))) user_io
(
	.clk_sys(clk_sys),
	.clk_sd(clk_sys),
	.conf_str(CONF_STR),

	.SPI_CLK(SPI_SCK),
	.SPI_SS_IO(CONF_DATA0),
	.SPI_MOSI(SPI_DI),
	.SPI_MISO(SPI_DO),

	.img_mounted(img_mounted),
	.img_size(img_size),
	.sd_conf(1'b0),
	.sd_ack_conf(),
	.sd_sdhc(1'b1),
	.sd_lba(sd_lba),
	.sd_rd(sd_rd),
	.sd_wr(sd_wr),
	.sd_ack(sd_ack),
	.sd_buff_addr(sd_buff_addr),
	.sd_din(sd_buff_din),
	.sd_dout(sd_buff_dout),
	.sd_dout_strobe(sd_buff_wr),

	.key_strobe(key_strobe),
	.key_code(key_code),
	.key_pressed(key_pressed),
	.key_extended(key_extended),

	.mouse_x(mouse_x),
	.mouse_y(mouse_y),
	.mouse_flags(mouse_flags),
	.mouse_strobe(mouse_strobe),

	.joystick_0(joystick_0),
	.joystick_1(joystick_1),

`ifdef USE_HDMI
	.i2c_start      (i2c_start      ),
	.i2c_read       (i2c_read       ),
	.i2c_addr       (i2c_addr       ),
	.i2c_subaddr    (i2c_subaddr    ),
	.i2c_dout       (i2c_dout       ),
	.i2c_din        (i2c_din        ),
	.i2c_ack        (i2c_ack        ),
	.i2c_end        (i2c_end        ),
`endif

	.buttons(buttons),
	.status(status),
	.scandoubler_disable(scandoubler_disable),
	.ypbpr(ypbpr),
	.no_csync(no_csync)

);
wire        ioctl_wr;
wire [24:0] ioctl_addr;
wire  [7:0] ioctl_dout;
wire        ioctl_download;
wire  [5:0] ioctl_index;
wire  [1:0] ioctl_ext_index;

data_io #(.START_ADDR(27'h080000)) data_io
(
	.clk_sys(clk_sys),

	.SPI_SCK(SPI_SCK),
	.SPI_SS2(SPI_SS2),
	.SPI_DI(SPI_DI),
	.SPI_DO(SPI_DO),

	.clkref_n(1'b0),
	.ioctl_wr(ioctl_wr),
	.ioctl_addr(ioctl_addr),
	.ioctl_dout(ioctl_dout),
	.ioctl_download(ioctl_download),
	.ioctl_index({ioctl_ext_index, ioctl_index})
);

///////////////////   CPU   ///////////////////
wire [15:0] addr;
wire  [7:0] cpu_din;
wire  [7:0] cpu_dout;
wire        nM1;
wire        nMREQ;
wire        nIORQ;
wire        nRD;
wire        nWR;
wire        nRFSH;
reg         nINT;
reg         reset;
wire        cold_reset = (mod[1] & Fn[11]) | init_reset;
wire        warm_reset =  mod[2] & Fn[11];
wire        port_we    = ~nIORQ & ~nWR & nM1;
wire        port_rd    = ~nIORQ & ~nRD & nM1;

always @(posedge clk_sys) reset <= buttons[1] | status[0] | cold_reset | warm_reset;
always @(posedge clk_sys) if(ce_cpu_p) nINT <= ~(INT_line | INT_frame | INT_midi);

T80pa cpu
(
	.RESET_n(~reset),
	.CLK(clk_sys),
	.CEN_p(ce_cpu_p),
	.CEN_n(ce_cpu_n),
	.WAIT_n(1),
	.INT_n(nINT),
	.NMI_n((mod[2:1]==0) & Fn[11]),
	.BUSRQ_n(1),
	.M1_n(nM1),
	.MREQ_n(nMREQ),
	.IORQ_n(nIORQ),
	.RD_n(nRD),
	.WR_n(nWR),
	.RFSH_n(nRFSH),
	.HALT_n(1),
	.A(addr),
	.DO(cpu_dout),
	.DI(cpu_din)
);

always_comb begin
	case({nMREQ, ~nM1 | nIORQ | nRD})
	    'b01: cpu_din = ext_ena ? ram_dout : 8'hFF;
	    'b10: cpu_din = asic_dout;
	 default: cpu_din = 8'hFF;
	endcase
end

reg init_reset = 1;
always @(posedge clk_sys) begin
	reg old_download;
	old_download <= ioctl_download;
	if(~ioctl_download & old_download & !ioctl_index) init_reset <= 0;
end


//////////////////   MEMORY   //////////////////
reg  [24:0] ram_addr;
always_comb begin
	casex({rom0_sel | rom1_sel, ext_ram, addr[15:14]})
		'b1X_XX: ram_addr = {5'h10,addr[15], addr[13:0]};
		'b01_X0: ram_addr = {ext_c_off,      addr[13:0]};
		'b01_X1: ram_addr = {ext_d_off,      addr[13:0]};
		'b00_00: ram_addr = {page_ab,        addr[13:0]};
		'b00_01: ram_addr = {page_ab + 1'b1, addr[13:0]};
		'b00_10: ram_addr = {page_cd,        addr[13:0]};
		'b00_11: ram_addr = {page_cd + 1'b1, addr[13:0]};
	endcase
end

wire       ram_busy;
wire [7:0] ram_dout;
sdram ram
(
	.*,
	.init(~locked),
	.clk(clk_mem),
	.addr(ram_addr),
	.dout(ram_dout),
	.din(cpu_dout),
	.we(~(rom0_sel | rom1_sel | ram_wp) & ~nMREQ & ~nWR & ext_ena),
	.rd(~nMREQ & ~nRD),

	.vid_addr1(vram_addr1),
	.vid_addr2(vram_addr2),
	.vid_data1(vram_dout1),
	.vid_data2(vram_dout2),

	.misc_addr(ioctl_addr),
	.misc_din(ioctl_dout),
	.misc_dout(),
	.misc_rd(0),
	.misc_we(ioctl_wr),
	.misc_busy()
);


//////////////////////  EXT.RAM  ////////////////////
reg  [7:0] ext_c;
reg  [7:0] ext_d;
wire [8:0] ext_c_off = 9'h40 + ext_c;
wire [8:0] ext_d_off = 9'h40 + ext_d;
reg        ext_dis;
wire       ext_ena   = ~(ext_ram & ext_dis);

wire       ext_c_sel = (addr[7:0] == 128);
wire       ext_d_sel = (addr[7:0] == 129);

always @(posedge clk_sys) begin
	reg old_we;
	old_we <= port_we;
	if(port_we & ~old_we) begin
		if(ext_c_sel) ext_c <= cpu_dout;
		if(ext_d_sel) ext_d <= cpu_dout;
	end

	if(reset) ext_dis <= status[5];
end


////////////////////  ASIC PORTS  ///////////////////
reg  [7:0] brdr;
wire [3:0] border_color = {brdr[5], brdr[2:0]};
wire       ear_out  = brdr[4];
wire       mic_out  = 0; // brdr[3]; it seems not used in SAM Coupe.

reg  [7:0] lmpr;
wire [4:0] page_ab  = lmpr[4:0];
wire       rom0_sel =~lmpr[5] & !addr[15:14];
wire       rom1_sel = lmpr[6] & &addr[15:14];
wire       ram_wp   = lmpr[7] & !addr[15:14];

reg  [7:0] hmpr;
wire [4:0] page_cd  = hmpr[4:0];
wire [1:0] mode3_hi = hmpr[6:5];
wire       ext_ram  = hmpr[7] &  addr[15];

wire       sid_sel  = (addr[7:0] == 212); // D4
wire       fdd_sel  = &addr[7:5] & ~addr[3]; // 224-231(E0-E7), 240-247(F0-F7)
wire       lptd_sel = (addr[7:0] == 232); // E8
wire       lpts_sel = (addr[7:0] == 233); // E9
//clut, hpen, lpen  = (addr[7:0] == 248); // F8
wire       stat_sel = (addr[7:0] == 249); // F9
wire       lmpr_sel = (addr[7:0] == 250); // FA
wire       hmpr_sel = (addr[7:0] == 251); // FB
//         vmpr_sel = (addr[7:0] == 252); // FC
wire       midi_sel = (addr[7:0] == 253); // FD
wire       kbdr_sel = (addr[7:0] == 254); // FE
wire       brdr_sel = (addr[7:0] == 254); // FE
//         attr_sel = (addr[7:0] == 255); // FF

always @(posedge clk_sys) begin
	reg old_we;
	
	if(reset) begin
		lmpr <= 0;
		hmpr <= 0;
		brdr <= 0;
	end else begin
		old_we <= port_we;
		if(port_we & ~old_we) begin
			if(brdr_sel) brdr <= cpu_dout;
			if(lmpr_sel) lmpr <= cpu_dout;
			if(hmpr_sel) hmpr <= cpu_dout;
		end
	end
end

reg [7:0] asic_dout;
always_comb begin
	casex({kbdr_sel, stat_sel, lmpr_sel, hmpr_sel, vid_sel, fdd_sel, kjoy_sel, lptd_sel | lpts_sel})
		'b1XXXXXXX: asic_dout = {soff, ~tape_in, 1'b0, hid_data};
		'b01XXXXXX: asic_dout = {key_data[7:5], ~INT_midi, ~INT_frame, 2'b11, ~INT_line};
		'b001XXXXX: asic_dout = lmpr;
		'b0001XXXX: asic_dout = hmpr;
		'b00001XXX: asic_dout = vid_dout;
		'b000001XX: asic_dout = fdd1_io ? fdd1_dout : fdd2_dout;
		'b0000001X: asic_dout = {2'b00, joystick_0[5:0] | joystick_1[5:0]};
		'b00000001: asic_dout = 0; // fake LPT port.
		'b00000000: asic_dout = 8'hFF;
	endcase
end

reg        INT_midi;
reg        midi_tx;
always @(posedge clk_sys) begin
	reg       old_we;
	reg [8:0] tx_time;
	reg       available;
	
	if(reset) begin
		INT_midi <=0;
		tx_time <= 0;
		midi_tx <= 0;
		available <= 0;
	end else begin
		if(ce_1m) begin
			if(tx_time) begin
				tx_time <= tx_time - 1'd1;
				if(tx_time == 15) INT_midi <= 1;
			end else begin
				{INT_midi, midi_tx} <= 0;
				if(available) begin
					midi_tx   <= 1;
					tx_time   <= 319;
					available <= 0;
				end
			end
		end

		old_we <= port_we;
		if(port_we & ~old_we & midi_sel) available <= 1;
	end
end

////////////////////   AUDIO   ///////////////////
wire [7:0] psg_ch_l;
wire [7:0] psg_ch_r;
wire       tape_in = 0; //tape is not implemented (yet?)

saa1099 psg
(
	.clk_sys(clk_sys),
	.ce(ce_psg),
	.rst_n(~reset),
	.cs_n((addr[7:0] != 255) | nIORQ),
	.a0(addr[8]),
	.wr_n(nWR),
	.din(cpu_dout),
	.out_l(psg_ch_l),
	.out_r(psg_ch_r)
);

wire [18:0] audioR = {1'b0, psg_ch_l, psg_ch_l, 2'd0} + {1'b0, ear_out, mic_out, tape_in, 15'd0} + {vox_l,vox_l,2'd0} + sid_out;
wire [18:0] audioL = {1'b0, psg_ch_r, psg_ch_r, 2'd0} + {1'b0, ear_out, mic_out, tape_in, 15'd0} + {vox_r,vox_r,2'd0} + sid_out;

sigma_delta_dac #(18) dac_l
(
	.CLK(clk_sys),
	.RESET(0),
	.DACin(audioL),
	.DACout(AUDIO_L)
);

sigma_delta_dac #(18) dac_r
(
	.CLK(clk_sys),
	.RESET(0),
	.DACin(audioR),
	.DACout(AUDIO_R)
);

`ifdef I2S_AUDIO
i2s i2s (
	.reset(1'b0),
	.clk(clk_sys),
	.clk_rate(32'd48_000_000),

	.sclk(I2S_BCK),
	.lrclk(I2S_LRCK),
	.sdata(I2S_DATA),

	.left_chan({~audioL[18], audioL[17:3]}),
	.right_chan({~audioR[18], audioR[17:3]})
);
`endif

`ifdef SPDIF_AUDIO
spdif spdif
(
	.clk_i(clk_sys),
	.rst_i(1'b0),
	.clk_rate_i(32'd48_000_000),
	.spdif_o(SPDIF),
	.sample_i({~audioR[18], audioR[17:3], ~audioL[18], audioL[17:3]})
);
`endif

reg [7:0] vox_l, vox_r;
always @(posedge clk_sys) begin
	reg old_we, old_stb;
	reg [7:0] data;

	if(reset) {vox_l, vox_r, old_stb, data} <= 0;
	else begin
		old_we <= port_we;
		if(port_we & ~old_we) begin
			if(lptd_sel) data <= cpu_dout;
			if(lpts_sel) begin
				if(~old_stb & cpu_dout[0]) vox_l <= data;
				if(old_stb & ~cpu_dout[0]) vox_r <= data;
				old_stb <= cpu_dout[0];
			end
		end
	end
end

// SID uses signed samples and requires special handling
softmuter sid_muter
(
	.clk_sys(clk_sys),
	.ce(ce_1m),

	.enable(sid_act),
	.vol_in({~sid_outraw[17], sid_outraw[16:0]}),
	.vol_out(sid_out)
);

wire [17:0] sid_out;
wire [17:0] sid_outraw;
wire        sid_act;

sid_top sid
(
	.clock(clk_sys),
	.reset(reset),
	.addr(addr[12:8]),
	.wren(port_we && sid_sel && video_mode), // disable in ZX mode
	.wdata(cpu_dout),

	.comb_wave_l(0),
	.comb_wave_r(0),
	.extfilter_en(1),

	.active(sid_act),
	.start_iter(ce_1m),
	.sample_left(sid_outraw)
);

////////////////////   VIDEO   ///////////////////
wire [18:0] vram_addr1;
wire [18:0] vram_addr2;
wire        vram_rd1;
wire        vram_rd2;
wire [15:0] vram_dout1;
wire [15:0] vram_dout2;
wire  [7:0] vid_dout;
wire        vid_sel;
wire        soff = (brdr[7] & video_mode[1]) | turbo_boot;
wire        INT_line;
wire        INT_frame;
wire        mem_contention;
wire        io_contention;
wire  [1:0] video_mode;
wire  [1:0] R,G,B;
wire        I;
wire        HBlank, VBlank, HSync, VSync;

video video
(
	.*,
	.wide(0),
	.full_zx(status[12:11] == 1),
	.din(cpu_dout),
	.dout(vid_dout),
	.dout_en(vid_sel)
);

mist_video #(.COLOR_DEPTH(4), .SD_HCNT_WIDTH(10), .OUT_COLOR_DEPTH(VGA_BITS), .BIG_OSD(BIG_OSD), .USE_BLANKS(1'b1)) mist_video (
	.clk_sys     ( clk_sys    ),

	// OSD SPI interface
	.SPI_SCK     ( SPI_SCK    ),
	.SPI_SS3     ( SPI_SS3    ),
	.SPI_DI      ( SPI_DI     ),

	// scanlines (00-none 01-25% 10-50% 11-75%)
	.scanlines   ( st_scanlines  ),

	// non-scandoubled pixel clock divider 0 - clk_sys/4, 1 - clk_sys/2
	.ce_divider  ( 3'd0       ),

	// 0 = HVSync 31KHz, 1 = CSync 15KHz
	.scandoubler_disable ( scandoubler_disable ),
	// disable csync without scandoubler
	.no_csync    ( no_csync   ),
	// YPbPr always uses composite sync
	.ypbpr       ( ypbpr      ),
	// Rotate OSD [0] - rotate [1] - left or right
	.rotate      ( 2'b00      ),
	// composite-like blending
	.blend       ( 1'b0       ),

	// video in
	.R           ( {R, R[1], I} ),
	.G           ( {G, G[1], I} ),
	.B           ( {B, B[1], I} ),

	.HSync       ( HSync      ),
	.VSync       ( VSync      ),
	.HBlank      ( HBlank     ),
	.VBlank      ( VBlank     ),

	// MiST video output signals
	.VGA_R       ( VGA_R      ),
	.VGA_G       ( VGA_G      ),
	.VGA_B       ( VGA_B      ),
	.VGA_VS      ( VGA_VS     ),
	.VGA_HS      ( VGA_HS     )
);

`ifdef USE_HDMI
i2c_master #(48_000_000) i2c_master (
	.CLK         (clk_sys),
	.I2C_START   (i2c_start),
	.I2C_READ    (i2c_read),
	.I2C_ADDR    (i2c_addr),
	.I2C_SUBADDR (i2c_subaddr),
	.I2C_WDATA   (i2c_dout),
	.I2C_RDATA   (i2c_din),
	.I2C_END     (i2c_end),
	.I2C_ACK     (i2c_ack),

	//I2C bus
	.I2C_SCL     (HDMI_SCL),
	.I2C_SDA     (HDMI_SDA)
);

mist_video #(.COLOR_DEPTH(4), .SD_HCNT_WIDTH(10), .OUT_COLOR_DEPTH(8), .BIG_OSD(BIG_OSD), .USE_BLANKS(1'b1), .VIDEO_CLEANER(1'b1)) hdmi_video (
	.clk_sys     ( clk_sys    ),

	// OSD SPI interface
	.SPI_SCK     ( SPI_SCK    ),
	.SPI_SS3     ( SPI_SS3    ),
	.SPI_DI      ( SPI_DI     ),

	// scanlines (00-none 01-25% 10-50% 11-75%)
	.scanlines   ( st_scanlines  ),

	// non-scandoubled pixel clock divider 0 - clk_sys/4, 1 - clk_sys/2
	.ce_divider  ( 3'd0       ),

	// 0 = HVSync 31KHz, 1 = CSync 15KHz
	.scandoubler_disable ( 1'b0 ),
	// disable csync without scandoubler
	.no_csync    ( 1'b0       ),
	// YPbPr always uses composite sync
	.ypbpr       ( 1'b0       ),
	// Rotate OSD [0] - rotate [1] - left or right
	.rotate      ( 2'b00      ),
	// composite-like blending
	.blend       ( 1'b0       ),

	// video in
	.R           ( {R, R[1], I} ),
	.G           ( {G, G[1], I} ),
	.B           ( {B, B[1], I} ),
	.HSync       ( HSync      ),
	.VSync       ( VSync      ),
	.HBlank      ( HBlank     ),
	.VBlank      ( VBlank     ),

	// MiST video output signals
	.VGA_R       ( HDMI_R      ),
	.VGA_G       ( HDMI_G      ),
	.VGA_B       ( HDMI_B      ),
	.VGA_VS      ( HDMI_VS     ),
	.VGA_HS      ( HDMI_HS     ),
	.VGA_DE      ( HDMI_DE     )
);

assign HDMI_PCLK = clk_sys;

`endif

////////////////////   HID   /////////////////////
wire [11:1] Fn;
wire  [2:0] mod;
wire  [7:0] key_data;
reg         autostart;
keyboard kbd( .*, .restart(rom0_sel & (addr == 0) & ~nMREQ & ~nRD));

wire        kjoy_sel = (addr[7:0] == 'h1F);
wire  [4:0] hid_data = key_data[4:0] & mouse_data
	& (addr[12] ? 5'b11111 : ~{joystick_0[1],  joystick_0[0], joystick_0[2], joystick_0[3], joystick_0[4] | joystick_0[5]})
	& (addr[11] ? 5'b11111 : ~{joystick_1[4] | joystick_1[5], joystick_1[3], joystick_1[2], joystick_1[0],  joystick_1[1]});

wire  [4:0] mouse_data;
mouse mouse( .*, .dout(mouse_data), .rd(kbdr_sel & &addr[15:8] & nM1 & ~nIORQ & ~nRD));


///////////////////   FDC   ///////////////////

reg fdd_num = 0;
always @(posedge clk_sys) begin
	if(sd_rd[1]|sd_wr[1]) fdd_num <= 1;
	if(sd_rd[0]|sd_wr[0]) fdd_num <= 0;
end

assign sd_buff_din = fdd_num ? fdd2_buf_dout : fdd1_buf_dout;
assign sd_lba      = fdd_num ? fdd2_lba      : fdd1_lba;



// FDD1
wire        fdd1_busy;
reg         fdd1_ready;
reg         fdd1_side;
wire        fdd1_io   = fdd_sel & ~addr[4] & ~nIORQ & nM1;
wire  [7:0] fdd1_dout;
wire  [7:0] fdd1_buf_dout;
wire [31:0] fdd1_lba;

always @(posedge clk_sys) begin
	reg old_wr;
	reg old_mounted;

	old_wr <= nWR;
	if(old_wr & ~nWR & fdd1_io) fdd1_side <= addr[2];

	old_mounted <= img_mounted[0];
	if(cold_reset) fdd1_ready <= 0;
		else if(~old_mounted & img_mounted[0]) fdd1_ready <= 1;
end

wd1793 #(1) fdd1
(
	.clk_sys(clk_sys),
	.ce(cpu_n),
	.reset(reset),
	.io_en(fdd1_io & fdd1_ready),
	.rd(~nRD),
	.wr(~nWR),
	.addr(addr[1:0]),
	.din(cpu_dout),
	.dout(fdd1_dout),

	.img_mounted(img_mounted[0]),
	.img_size(img_size[19:0]),
	.sd_lba(fdd1_lba),
	.sd_rd(sd_rd[0]),
	.sd_wr(sd_wr[0]),
	.sd_ack(sd_ack),
	.sd_buff_addr(sd_buff_addr),
	.sd_buff_dout(sd_buff_dout),
	.sd_buff_din(fdd1_buf_dout),
	.sd_buff_wr(sd_buff_wr),

	.wp(~status[4]),

	.size_code(4),
	.layout(ioctl_ext_index == 2),
	.side(fdd1_side),
	.ready(fdd1_ready),
	.prepare(fdd1_busy),

	.input_active(0),
	.input_addr(0),
	.input_data(0),
	.input_wr(0),
	.buff_din(0)
);

always @(posedge clk_sys) begin
	reg old_busy;
	integer counter;

	if(reset) begin
		autostart <= 0;
		counter   <= 0;
	end else begin
		old_busy <= fdd1_busy;
		if(old_busy & ~fdd1_busy & fdd1_ready) begin
			autostart <= 1;
			counter   <= 1000000;
		end else if(ce_6mp && counter) begin
			counter   <= counter - 1;
			if(counter == 1) autostart <= 0;
		end
	end
end



// FDD2
reg         fdd2_ready;
reg         fdd2_side;
wire        fdd2_io   = fdd_sel & addr[4] & ~nIORQ & nM1;
wire  [7:0] fdd2_dout;
wire  [7:0] fdd2_buf_dout;
wire [31:0] fdd2_lba;

always @(posedge clk_sys) begin
	reg old_wr;
	reg old_mounted;

	old_wr <= nWR;
	if(old_wr & ~nWR & fdd2_io) fdd2_side <= addr[2];

	old_mounted <= img_mounted[1];
	if(cold_reset) fdd2_ready <= 0;
		else if(~old_mounted & img_mounted[1]) fdd2_ready <= 1;
end

wd1793 #(1) fdd2
(
	.clk_sys(clk_sys),
	.ce(cpu_n),
	.reset(reset),
	.io_en(fdd2_io & fdd2_ready),
	.rd(~nRD),
	.wr(~nWR),
	.addr(addr[1:0]),
	.din(cpu_dout),
	.dout(fdd2_dout),

	.img_mounted(img_mounted[1]),
	.img_size(img_size[19:0]),
	.sd_lba(fdd2_lba),
	.sd_rd(sd_rd[1]),
	.sd_wr(sd_wr[1]),
	.sd_ack(sd_ack),
	.sd_buff_addr(sd_buff_addr),
	.sd_buff_dout(sd_buff_dout),
	.sd_buff_din(fdd2_buf_dout),
	.sd_buff_wr(sd_buff_wr),

	.wp(~status[4]),

	.size_code(4),
	.layout(ioctl_ext_index == 2),
	.side(fdd2_side),
	.ready(fdd2_ready),

	.input_active(0),
	.input_addr(0),
	.input_data(0),
	.input_wr(0),
	.buff_din(0)
);

endmodule

module softmuter
(
	input         clk_sys,
	input         ce,
	
	input         enable,
	input  [17:0] vol_in,
	output [17:0] vol_out
);

reg [17:0] att = '1;
assign vol_out = (vol_in > att) ? vol_in - att : 18'd0;

always @(posedge clk_sys) begin
	if(ce) begin
		if( enable &&   att) att <= att - 1'd1;
		if(~enable && ~&att) att <= att + 1'd1;
	end
end

endmodule
