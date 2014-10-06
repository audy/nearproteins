#!/usr/bin/env ruby

require 'dna'
require 'socket'
require 'json'

HOSTNAME = '127.0.0.1'
PORT = '1234'

socket = TCPSocket.open HOSTNAME, PORT

File.open('proteins.fasta') do |handle|
  records = Dna.new handle, :format => :fasta

  records.each do |record|
    socket.puts "GET #{record.sequence}"
    resp = JSON.parse(socket.gets)
    p resp
  end
end
