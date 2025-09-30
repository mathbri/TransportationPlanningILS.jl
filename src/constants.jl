# This file contains all constants used in the project

# Actual values are Symbol, those list are here for verification
const NODE_TYPES = [:supplier, :plant, :platform, :pol, :pod]
const COMMON_NODE_TYPES = [:platform, :pol, :pod]
const ARC_TYPES = [:direct, :outsource, :cross_plat, :delivery, :oversea, :shortcut]
const BP_ARC_TYPES = [:outsource, :direct, :cross_plat, :delivery, :oversea]
const COMMON_ARC_TYPES = [:cross_plat, :delivery, :oversea]

const EPS = 1e-5
const INFINITY = 1e9
const VOLUME_FACTOR = 100

const LAND_CAPACITY = 7000
const SEA_CAPACITY = 6500

const PERTURBATIONS = [:single_plant, :attract_reduce]
const MAX_MILP_VAR = 2_000_000
const MILP_TIME_LIMIT = 120

const NETWORK_DESIGN_COLUMNS = [
    "route_id",
    "supplier_account",
    "customer_account",
    "part_number",
    "packaging",
    "quantity_part_in_route",
    "delivery_date",
    "point_account",
    "point_index",
    "point_date",
    "shipment_id",
]

const SHIPMENT_INFO_COLUMNS = [
    "shipment_id",
    "source_point_account",
    "destination_point_account",
    "point_start_date",
    "point_end_date",
    "type",
    "volume",
    "transport_cost",
    "carbon_cost",
    "platform_cost",
]

const SHIPMENT_CONTENT_COLUMNS = [
    "content_id",
    "shipment_id",
    "part_number",
    "part_supplier_account",
    "part_customer_account",
    "quantity",
    "packaging_size",
    "volume",
]

const ROUTE_COLUMNS = [
    "route_id",
    "supplier_account",
    "customer_account",
    "point_account",
    "point_number",
    "point_type",
]